{-# language FlexibleContexts    #-}
{-# language RebindableSyntax    #-}
{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}

module Gridding where

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A
import Data.Array.Accelerate.Math.DFT.Centre              as A
import Data.Array.Accelerate.Math.FFT                     as A

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..), maybe)

import Control.Lens as L (_1, _2, _3, _4, _5)
import Data.Maybe (fromJust, fromMaybe)
--import Data.Function.Memoize
import System.IO.Unsafe (unsafePerformIO)
import Debug.Trace

data Kwargs = Kwargs { patHorShift :: Maybe Int
                     , patVerShift :: Maybe Int
                     , patTransMat :: Maybe (Acc (Matrix Double))
                     , wstep :: Maybe Int
                     , qpx :: Maybe Int              -- The oversampling of the convolution kernel
                     , npixFF :: Maybe Int
                     , npixKern :: Maybe Int         -- Kernel size (if it is same size in x and y)
                     }

data OtherImagingArgs = OtherImagingArgs
                      { convolutionKernel :: Maybe ((Array DIM4 (Complex Double)))
                      , kernelCache :: Maybe KernelF
                      , kernelFunction :: Maybe KernelF
                      }

type KernelF = Exp Double                            -- Field of view
             -> Acc (Scalar Double)                  -- Baseline distance to the projection plane
             -> Maybe (Exp Int)                      -- The id of the antenna 1
             -> Maybe (Exp Int)                      -- The id of the antenna 2
             -> Maybe (Exp Double)                   -- The time
             -> Maybe (Exp Double)                   -- The frequency
             -> Kwargs                               -- Additional kernel arguments
             -> Acc (Array DIM4 (Complex Double))    -- [Qpx,Qpx,s,s] shaped oversampled convolution kernels

type WKernelF = Exp Double                            -- Field of view
              -> Acc (Scalar Double)                  -- Baseline distance to the projection plane
              -> Kwargs                               -- Additional kernel arguments
              -> Acc (Array DIM4 (Complex Double))    -- [Qpx,Qpx,s,s] shaped oversampled convolution kernels

type AKernelF = Exp Double                            -- Field of view
              -> Exp Int                              -- The id of the antenna
              -> Exp Double                           -- The time
              -> Exp Double                           -- The frequency
              -> Acc (Matrix (Complex Double))        -- The a-kernel


noArgs :: Kwargs
noArgs = Kwargs Nothing Nothing Nothing Nothing Nothing Nothing Nothing

noOtherArgs :: OtherImagingArgs
noOtherArgs = OtherImagingArgs Nothing Nothing Nothing
-------------------------------
-- Gridding Accelerate code
type ImagingFunction = Double                                        -- (theta) Field of view size
                     -> Int                                          -- (lam) Grid size
                     -> Acc (Vector (Double, Double, Double))        -- (uvw) all the uvw baselines (coordinates) (lenght : n * (n-1))
                     -> Acc (Vector (Int, Int, Double, Double))      -- (src) (Antenna 1, Antenna 2, The time (in MJD UTC), Frequency (Hz)) (lenght : n * (n-1))
                     -> Acc (Vector (Complex Double))                -- (vis) visibility  (length n * (n-1))
                     -> Kwargs                                       -- (kwargs) The additional options for the kernel
                     -> OtherImagingArgs                             -- Some other arguments that not all imaging functions use
                     -> Acc (Matrix (Complex Double))

-- # Simple imaging
simple_imaging :: ImagingFunction
simple_imaging theta lam uvw src vis kwargs otargs =
    let lamf = fromIntegral lam
        n = P.round(theta * lamf)
        ne = constant n

        p = map (`div3` constant lamf) uvw
        czero = constant 0
        guv = fill (index2 ne ne) czero
    in grid guv p vis

grid :: Acc (Matrix (Complex Double))                    -- Destination grid N x N of complex numbers
     -> Acc (Vector (Double, Double, Double))            -- The uvw baselines, but between -.5 and .5
     -> Acc (Vector (Complex Double))                    -- The visiblities
     -> Acc (Matrix (Complex Double))
grid a p v = permute (+) a (\i -> xy ! i) v
    where
        Z :. n :. _ = unlift (shape a) :: Z :. Exp Int :. Exp Int
        halfn = n `div` 2
        nf = A.fromIntegral n :: Exp Double
        xy = map gridxy p :: Acc (Vector DIM2)

        -- Note the order, python specified as a[y, x] += v
        -- We do the same. This is correct, first y then x in the index2 function
        gridxy :: Exp (Double, Double, Double) -> Exp DIM2
        gridxy (unlift -> (x, y, _ :: Exp Double)) = index2 (toGridCell y) (toGridCell x)

        toGridCell :: Exp Double -> Exp Int
        toGridCell f = halfn + floor (0.5 + nf * f)

-- # Normal imaging
conv_imaging :: ImagingFunction
conv_imaging theta lam uvw src vis kwargs otargs =
    let kv = fromJust $ convolutionKernel otargs
        lamf = fromIntegral lam
        n = P.round(theta * lamf)
        ne = constant n

        p = map (`div3` constant lamf) uvw
        czero = constant 0
        guv = fill (index2 ne ne) czero
    in convgrid (use kv) guv p vis

frac_coord :: Exp Int                                   -- The height/width
           -> Exp Int                                   -- Oversampling factor
           -> Acc (Vector Double)                       -- a uvw baseline, but between -.5 and .5
           -> Acc (Vector (Int, Int))
frac_coord n qpx p =
    let halfn = n `div` 2
        halfnf = A.fromIntegral halfn
        nf = A.fromIntegral n
        qpxf = A.fromIntegral qpx
        qpxfrac = 0.5 / qpxf

        x   = map (\a -> halfnf + a * nf) p
        flx = map (\a -> floor (a + qpxfrac) ) x
        fracx = zipWith (\a b -> round ((a - (A.fromIntegral b) )* qpxf) ) x flx
    in zip flx fracx

frac_coords :: Exp (Int, Int)                           -- (Heigt, width) of the grid
            -> Exp Int                                  -- Oversampling factor
            -> Acc (Vector (Double, Double, Double))    -- The uvw baselines, but between -.5 and .5
            -> Acc (Vector (Int, Int, Int, Int))
-- NOTE we should give it in height, width order.
frac_coords (unlift -> (h, w) ) qpx p =
    let (u, v, _) = unzip3 p
        (x, xf)   = unzip (frac_coord w qpx u)
        (y, yf)   = unzip (frac_coord h qpx v)
    in zip4 x xf y yf

convgrid :: Acc (Array DIM4 (Complex Double))        -- The oversampled convolution kernel
         -> Acc (Matrix (Complex Double))            -- Destination grid N x N of complex numbers
         -> Acc (Vector (Double, Double, Double))    -- The uvw baselines, but between -.5 and .5
         -> Acc (Vector (Complex Double))            -- The visiblities
         -> Acc (Matrix (Complex Double))
convgrid gcf a p v =
    let Z :. qpx :. _ :. gh :. gw = unlift (shape gcf) :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int
        Z :. height :. width = unlift (shape a) :: Z :. Exp Int :. Exp Int
        halfgh = gh `div` 2
        halfgw = gw `div` 2
        cc0 = (lift (constant 0.0 :+ constant 0.0))

        coords = frac_coords (lift (height, width)) qpx p
        (cx, cxf, cy, cyf) = unzip4 coords

        -- Shift the x and y -0.5*gw and -0.5*gh
        -- From here on the convolution kernel maps to the destination grid
        newcx = map (\x -> x - halfgw) cx
        newcy = map (\y -> y - halfgh) cy
        coordsAndVis = zip5 newcx cxf newcy cyf v
        -- We replicate gh * gw times, because each visibility is added is added to so many points,
        -- centered around the orignal x and y
        coordsAndVisRep = replicate (lift (Z :. All :. gh  :. gw)) coordsAndVis

        -- Add the correct offset, and multiply vis with the kernel
        withkern = imap getComplexAndAddOffset coordsAndVisRep
        -- Fix values that are now set out of bounds
        fixedBounds = map fixer withkern
        (x,y, val) = unzip3 fixedBounds

        fixer = fixoutofbounds width height cc0

        indexer id =
            let y' = y ! id
                x' = x ! id
            in index2 y' x'

        getComplexAndAddOffset :: Exp DIM3 -> Exp (Int, Int, Int, Int, Complex Double) -> Exp (Int, Int, Complex Double)
        getComplexAndAddOffset
            (unlift . unindex3 -> (_, i, j) ::(Exp Int,Exp Int,Exp Int))
            (unlift -> (x, xf, y, yf, vis)::(Exp Int,Exp Int,Exp Int,Exp Int,Exp (Complex Double))) =
            lift ( x + j
                 , y + i
                 , vis * gcf ! lift (Z :. yf :. xf :. i :. j) )
    in permute (+) a indexer val

convgrid2 :: Acc (Array DIM5 (Complex Double))        -- The oversampled convolution kernel
          -> Acc (Matrix (Complex Double))            -- Destination grid N x N of complex numbers
          -> Acc (Vector (Double, Double, Double))    -- The uvw baselines, but between -.5 and .5
          -> Acc (Vector Int)                         -- *DIF* The wbin index of the convolution kernel
          -> Acc (Vector (Complex Double))            -- The visiblities
          -> Acc (Matrix (Complex Double))
convgrid2 gcf a p wbin v =
    let Z :. _ :. qpx :. _ :. gh :. gw = unlift (shape gcf) :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int :. Exp Int
        Z :. height :. width = unlift (shape a) :: Z :. Exp Int :. Exp Int
        halfgh = gh `div` 2
        halfgw = gw `div` 2
        cc0 = (lift (constant 0.0 :+ constant 0.0))

        coords = frac_coords (lift (height, width)) qpx p
        (cx, cxf, cy, cyf) = unzip4 coords

        -- Shift the x and y -0.5*gw and -0.5*gh
        -- From here on the convolution kernel maps to the destination grid
        newcx = map (\x -> x - halfgw) cx
        newcy = map (\y -> y - halfgh) cy
        coordsAndVis = zip6 wbin newcx cxf newcy cyf v
        -- We replicate gh * gw times, because each visibility is added is added to so many points,
        -- centered around the orignal x and y
        coordsAndVisRep = replicate (lift (Z :. All :. gh  :. gw)) coordsAndVis

        -- Add the correct offset, and multiply vis with the kernel
        withkern = imap getComplexAndAddOffset coordsAndVisRep
        -- Fix values that are now set out of bounds
        fixedBounds = map fixer withkern
        (x,y, val) = unzip3 fixedBounds

        fixer = fixoutofbounds width height cc0

        indexer id =
            let y' = y ! id
                x' = x ! id
            in index2 y' x'

        getComplexAndAddOffset :: Exp DIM3 -> Exp (Int, Int, Int, Int, Int, Complex Double) -> Exp (Int, Int, Complex Double)
        getComplexAndAddOffset
            (unlift . unindex3 -> (_, i, j) ::(Exp Int,Exp Int,Exp Int))
            (unlift -> (wbin, x, xf, y, yf, vis)::(Exp Int,Exp Int,Exp Int,Exp Int,Exp Int,Exp (Complex Double))) =
            lift ( x + j
                 , y + i
                 , vis * gcf ! lift (Z :. wbin :. yf :. xf :. i :. j) )
    in permute (+) a indexer val

-- Imaging with a w kernel TODO: It differs from the normal implementation atm, it will only work with w-kernels, not aw-kernels
w_cache_imaging :: ImagingFunction
w_cache_imaging theta lam uvw src vis
 kwargs@Kwargs{wstep = wstep', qpx = qpx'}
 otargs@OtherImagingArgs{kernelCache = kernel_cache', kernelFunction = kernel_fn'} =
    let
        -- They use a LRU cache for the function, to memoize the results. Maybe usefull. Not sure yet.
        -- So maybe use LRU cache implementation (lrucache package) or memoization (memoize package)
        -- For now we just call the function.
        default_kernelf theta w _ _ _ _ kwargs = w_kernel theta w kwargs
        kernel_fn    = fromMaybe default_kernelf kernel_fn'
        kernel_cache = fromMaybe kernel_fn kernel_cache'
        only_wcache theta w kwargs = kernel_cache theta w Nothing Nothing Nothing Nothing kwargs
        wstep_ = fromMaybe 2000 wstep'
        wstep = constant $ wstep_

        lamf = fromIntegral lam
        n = P.round(theta * lamf)
        ne = constant n

        p = map (`div3` constant lamf) uvw
        czero = constant 0
        guv = fill (index2 ne ne) czero :: Acc (Matrix (Complex Double))


        -- We just make them all here, that's easiest for now
        --thekernels = fromList (Z :. wstep) [0..] :: Vector Double
        (_,_,w) = unzip3 uvw
        roundedw = map (\w' -> wstep * round ( w' / A.fromIntegral wstep)) w
        maxw = the $ maximum roundedw
        minw = the $ minimum roundedw
        steps = ((maxw - minw) `div` wstep ) + 1
        (cpumin, cpumax, cpusteps) = CPU.run (unit $ lift (minw, maxw, steps)) `indexArray` Z

        wbins = map (\rw' -> (rw' - constant cpumin) `div` wstep) roundedw
        
        makeWKernel :: Int -> Acc (Array DIM5 (Complex Double))
        makeWKernel i = let wbin = fromList Z $ (fromIntegral $ i * wstep_ + cpumin) : []
                        in use $ compiledMakeWKernel wbin
        
        compiledMakeWKernel = CPU.runN makeWKernel'
        
        makeWKernel' :: Acc (Scalar Double) -> Acc (Array DIM5 (Complex Double))
        makeWKernel' wbin = make5D . map conjugate $ only_wcache (constant theta) wbin kwargs

        
        make5D mat = let (Z :. yf :. xf :. y :. x) =  (unlift . shape) mat :: (Z :.Exp Int :. Exp Int :. Exp Int :. Exp Int)
                         newsh = lift (Z :. constant 1 :. yf :. xf :. y :. x) :: Exp DIM5
                     in reshape newsh mat

        thekernels = CPU.runN $ compute $ myfor (cpusteps - 1) (\i old -> concatOn _5 (makeWKernel i) old) (makeWKernel (cpusteps - 1))
    in convgrid2 (use thekernels) guv p wbins vis

-- Imaging with aw-caches
aw_cache_imaging :: ImagingFunction
aw_cache_imaging theta lam uvw src vis
 kwargs@Kwargs{wstep = wstep', qpx = qpx'}
 otargs@OtherImagingArgs{kernelCache = kernel_cache'} =
    let
        -- They use a LRU cache for the function, to memoize the results. Maybe usefull. Not sure yet.
        -- So maybe use LRU cache implementation (lrucache package) or memoization (memoize package)
        -- For now we just call the function.
        kernel_cache = fromJust kernel_cache'
        only_wcache theta w kwargs = kernel_cache theta w Nothing Nothing Nothing Nothing kwargs
        wstep_ = fromMaybe 2000 wstep'
        wstep = constant $ wstep_

        lamf = fromIntegral lam
        n = P.round(theta * lamf)
        ne = constant n

        p = map (`div3` constant lamf) uvw
        czero = constant 0
        guv = fill (index2 ne ne) czero :: Acc (Matrix (Complex Double))


        -- We just make them all here, that's easiest for now
        --thekernels = fromList (Z :. wstep) [0..] :: Vector Double
        (_,_,w) = unzip3 uvw
        roundedw = map (\w' -> wstep * round ( w' / A.fromIntegral wstep)) w
        maxw = the $ maximum roundedw
        minw = the $ minimum roundedw
        l = length roundedw
        steps = ((maxw - minw) `div` wstep ) + 1
        (cpumin, cpumax, cpusteps, cpul) = CPU.run (unit $ lift (minw, maxw, steps, l)) `indexArray` Z

        wbins = enumFromN (shape w) c0
    
        makeAWKernel :: Int -> Acc (Array DIM5 (Complex Double))
        makeAWKernel i = 
            let i' = use $ fromList Z [i]
                id = index1 (the i')
                wbin = unit . A.fromIntegral $ roundedw ! id
                (a1, a2, t, f) = unlift (src ! id) :: (Exp Int, Exp Int, Exp Double, Exp Double)
            in compute $ make5D . map conjugate $ kernel_cache (constant theta) wbin (Just a1) (Just a2) (Just t) (Just f) kwargs
        make5D mat = let (Z :. yf :. xf :. y :. x) =  (unlift . shape) mat :: (Z :.Exp Int :. Exp Int :. Exp Int :. Exp Int)
                         newsh = lift (Z :. constant 1 :. yf :. xf :. y :. x) :: Exp DIM5
                     in reshape newsh mat

        thekernels = myfor (cpul - 1) (\i old -> concatOn _5 (makeAWKernel i) old) (makeAWKernel (cpul - 1))
    in convgrid2 thekernels guv p wbins vis
    
----------------------------------
-- Processing the imaging functions
do_imaging :: Double                                 -- Field of view size
           -> Int                                    -- Grid size
           -> Vector (Double, Double, Double)        -- all the uvw baselines (coordinates) (lenght : n * (n-1))
           -> Vector (Int, Int, Double, Double)      -- (Antenna 1, Antenna 2, The time (in MJD UTC), Frequency (Hz)) (lenght : n * (n-1))
           -> Vector (Complex Double)                -- visibility  (length n * (n-1))
           -> ImagingFunction
           -> Kwargs
           -> OtherImagingArgs
           -> Acc (Matrix Double,Matrix Double, Scalar Double)
do_imaging theta lam uvw src vis imgfn kwargs otargs =
    let
        -- TODO: if src == None: src = numpy.ndarray((len(vis), 0))
        len = arrayShape vis
        uvw0 = use uvw
        vis0 = use vis
        src0 = use src
        -- Mirror baselines such that v>0
        (uvw1,vis1) = unzip $ mirror_uvw uvw0 vis0

        -- Determine weights
        ones = fill (lift len) 1
        wt = doweight theta lam uvw1 ones

        -- Make image
        cdrt = imgfn theta lam uvw1 src0 (zipWith (*) wt vis1) kwargs otargs
        drt = (map real . ifft . make_grid_hermitian) cdrt
        -- Make point spread function
        c = imgfn theta lam uvw1 src0 wt kwargs otargs
        psf =  (map real . ifft . make_grid_hermitian) c
        -- Normalise
        pmaxv = (maximum . flatten) psf
        pmax = the pmaxv
        res0 = map (/ pmax) drt
        res1 = map (/ pmax) psf
        res2 = pmaxv
    in lift (res0, res1, res2)

mirror_uvw :: Acc (Vector (Double, Double, Double))                     -- all the uvw baselines (coordinates) (lenght : n * (n-1))
           -> Acc (Vector (Complex Double))                             -- visibility  (length n * (n-1))
           -> Acc (Vector ((Double, Double, Double), Complex Double))
mirror_uvw uvw vis =
    let
        (u,v,w) = unzip3 uvw
        uvwvis = zip4 u v w vis
        mirrorv o@(unlift -> (u,v,w,vis) :: (Exp Double, Exp Double, Exp Double, Exp (Complex Double))) =
            if v < 0
                then lift ((-u, -v, -w), conjugate vis)
                else lift ((u,v,w),vis)
    in map mirrorv uvwvis

doweight :: Double                                                  -- Field of view size
         -> Int                                                     -- Grid size
         -> Acc (Vector (Double, Double, Double))                   -- all the uvw baselines (coordinates) (lenght : n * (n-1))
         -> Acc (Vector (Complex Double))                           -- Array of ones (same size as visibility  (length n * (n-1)))
         -> Acc (Vector (Complex Double))
doweight theta lam p v =
    let
        n = P.round(theta * lamf)
        lamf = fromIntegral lam
        ne = constant n :: Exp Int
        lamfe = constant lamf
        gw = fill (index2 ne ne) 0 :: Acc (Matrix Double)
        coords = frac_coords (lift (ne, ne)) 1 (map (`div3` lamfe) p)
        (x,_,y,_) = unzip4 coords
        ones  = fill (shape x) 1 :: Acc (Vector Double)
        xyindex = zipWith index2 y x
        weights = permute (+) gw (\id -> xyindex ! id) ones

        newv = imap (\id v -> v / lift (weights ! (xyindex ! id) :+ 0)) v
    in newv

make_grid_hermitian :: Acc (Matrix (Complex Double)) -> Acc (Matrix (Complex Double))
make_grid_hermitian guv = let
        sh = shape guv
        Z :. n :. _ = unlift sh :: Z :. Exp Int :. Exp Int
        indexer (unlift -> Z :. y :. x) = if x== 0 || y == 0
            then index2 y x
            else index2 (n - y) (n - x)
        mapper (unlift -> Z :. y :. x) v = if x== 0 || y == 0
            then 0
            else conjugate v

        -- Make mirror image, then add its conjugate to the original grid.
        -- We mirror on the zero point, which is off-center if the grid has even size
        oddReshape = (reverseOn _2 . reverseOn _1) guv
        evenReshape = backpermute sh indexer guv

        oddConj = map conjugate oddReshape
        evenConj = imap mapper evenReshape

        reshapedConj = if even n then evenConj else oddConj
    in zipWith (+) guv reshapedConj


----------------------
-- Kernels
w_kernel :: WKernelF
w_kernel theta w
    kwargs@Kwargs{npixFF = npixFF', npixKern = npixKern', qpx = qpx'} =
    let
        npixFF = constant $ fromJust npixFF'
        npixKern = constant $ fromJust npixKern'
        qpx = constant $ fromJust qpx'
        (l, m) = unlift $ kernel_coordinates npixFF theta kwargs :: (Acc (Matrix Double), Acc (Matrix Double))
        kern = w_kernel_function l m w
    in kernel_oversample kern npixFF qpx npixKern

kernel_coordinates :: Exp Int -> Exp Double -> Kwargs -> Acc ((Matrix Double), (Matrix Double))
kernel_coordinates n theta kwargs =
    let
        dl = (constant . fromIntegral . fromMaybe 0 . patHorShift) kwargs
        dm = (constant . fromIntegral . fromMaybe 0 . patVerShift) kwargs
        (l,m) = unlift $ coordinates2 n :: (Acc (Matrix Double), Acc (Matrix Double))
        lm1 = map (liftTupf (*theta) (*theta)) (zip l m)
        lm2 = case patTransMat kwargs of
            Nothing -> lm1
            Just t  -> map (\pair -> let (x, y) = unlift pair in lift
                ( (t ! index2 0 0) * x + (t ! index2 1 0) * y
                , (t ! index2 0 1) * x + (t ! index2 1 1) * y ) ) lm1

        lm3 = map (liftTupf (+dl) (+dm)) lm2
    in lift $ unzip lm3

coordinates2 :: Exp Int -> Acc (Matrix Double, Matrix Double)
coordinates2 n =
    let
        n2 = n `div` 2
        sh = index1 n
        sh2 = index2 n n
        start = A.fromIntegral (-n2) * step
        step = 1 / A.fromIntegral n
        base = enumFromStepN sh start step :: Acc (Vector Double)
        samecolumns = replicate (lift (Z :. n :. All)) base
        samerows = replicate (lift (Z :. All :. n)) base
    in lift (samecolumns, samerows)


w_kernel_function :: Acc (Matrix Double)              -- Horizontal image coordinates
                  -> Acc (Matrix Double)              -- Vertical image coordinates
                  -> Acc (Scalar Double)              -- Baseline distance to the projection plane
                  -> Acc (Matrix (Complex Double))    -- N x N array with the far field
w_kernel_function l m w' =
    let r2 = zipWith (\x y -> x*x + y*y) l m
    {-  dl and dm seem to be always zero anyway in the reference code
        dl = 0
        dm = 0
        ph1 = zipWith3 (\r2' l' m' -> 1 - sqrt (1 - r2') - dl * l' - dm * m') r2 l m
        ph2 = map (\r2' -> 1 - sqrt (1 - r2')) r2
        ph  = if dl == 0 && dm == 0 then ph2 else ph1
    -}
        w  = the w'
        ph = map (\r2' -> 1 - sqrt (1 - r2')) r2
        cp = map (\ph' -> (exp . lift) (0 :+ 2 * pi * w * ph') ) ph
    in cp

kernel_oversample :: Acc (Matrix (Complex Double))       -- Far field pattern
                  -> Exp Int                             -- Image size without oversampling
                  -> Exp Int                             -- Factor to oversample by -- there will be Qpx x Qpx convolution arl
                  -> Exp Int                             -- Size of convolution function to extract
                  -> Acc (Array DIM4 (Complex Double))   -- array of shape [ov, ou, v, u], e.g. with sub-pixel
kernel_oversample ff n qpx s =
    let
        -- Pad the far field to the required pixel size
        padff = pad_mid ff (n*qpx)
        -- Obtain oversampled uv-grid
        af = ifft padff
    in extract_oversampled af qpx s

pad_mid :: Acc (Matrix (Complex Double))    -- The input far field. Should be smaller than NxN
        -> Exp Int                          -- The desired far field size
        -> Acc (Matrix (Complex Double))    -- The far field, NxN
pad_mid ff n =
    let
        Z :. n0 :. n0w = (unlift . shape) ff :: Z :. Exp Int :. Exp Int
        result = if n == n0 then ff else padded
        pad_width = lift ((n `div` 2)-(n0 `div` 2), ((n+1) `div` 2)-((n0+1) `div` 2))
        padded = padder ff pad_width pad_width 0
    in result

extract_oversampled :: Acc (Matrix (Complex Double))        -- grid from which to extract
                    -> Exp Int                              -- oversampling factor
                    -> Exp Int                              -- size of section
                    -> Acc (Array DIM4 (Complex Double))    -- Return the oversampled kernels with shape qpx x qpx x n x n
extract_oversampled a qpx n =
    let
        -- In contrast to the python code, we return all the oversampled kernels at once, instead of one at a time
        na = (indexHead . shape) a
        cons = (na `div` 2) - qpx * (n `div` 2)
        m x = cons - x

        news = lift (Z :. qpx :. qpx :. n :. n) :: Exp DIM4
        indexer (unlift -> Z :. yf :. xf :. y :. x :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int)
            =  let newy = m yf + qpx * y
                   newx = m xf + qpx * x
               in index2 newy newx

        w_kern = backpermute news indexer a
        qpx2 = A.fromIntegral (qpx * qpx)
    in map (*qpx2) w_kern 

aw_kernel_fn :: AKernelF 
             -> Maybe WKernelF
             -> KernelF
aw_kernel_fn a_kernel_fn w_kernel_fn' theta w a1' a2' t' f'
    kwargs@Kwargs{qpx = qpx', npixKern = n'} =
    let 
        w_kernel_fn = fromMaybe w_kernel w_kernel_fn'
        qpx = fromJust qpx'
        a1 = fromJust a1'
        a2 = fromJust a2'
        t  = fromJust t'
        f  = fromJust f'
        n  = fromJust n'

        a1kern = a_kernel_fn theta a1 t f
        a2kern = a_kernel_fn theta a2 t f

        akern = convolve2d a1kern a2kern

        wkern = w_kernel_fn theta w kwargs

        getkern :: Elt a => Exp Int -> Exp Int -> Acc (Array DIM4 a) -> Acc (Matrix a)
        getkern yf xf a = slice a (lift (Z :. yf :. xf :. All :. All))

        newd = constant $ Z :. qpx :. qpx :. n :. n

        awkernf yf xf =  reshape newd $ convolve2d akern (getkern (constant yf) (constant xf) wkern)
        temp y = myfor (qpx - 1) (\x old -> concatOn _3 (awkernf y x) old) (awkernf y (qpx - 1))

        result = myfor (qpx - 1) (\y old -> concatOn _4 (temp y) old) (temp (qpx - 1))

    in result

-- the kernels for Aw-gridding will normally be chosen to /not/ overflow the borders 
-- Thus we chose the simpler method
convolve2d :: Acc (Matrix (Complex Double)) -> Acc (Matrix (Complex Double)) -> Acc (Matrix (Complex Double))
convolve2d a1 a2 = 
    let
        --Actually I expect 
        (Z :. n :. _) = (unlift . shape) a1 :: Z :. Exp Int :. Exp Int
        n2 = A.fromIntegral $ n * n

        a1fft = fft2D Inverse . ishift2D $ a1
        a2fft = fft2D Inverse . ishift2D $ a2
        a1a2 = zipWith (*) a1fft a2fft
        convolved = shift2D . fft2D Forward $ a1a2
    in map (*n2) convolved


----------------------
-- Fourier transformations
fft :: Acc (Matrix (Complex Double)) -> Acc (Matrix (Complex Double))
fft = shift2D . fft2D Forward . ishift2D

ifft :: Acc (Matrix (Complex Double)) -> Acc (Matrix (Complex Double))
ifft = shift2D . fft2D Inverse . ishift2D

------------------------
-- Helper functions
div3 :: Exp (Double, Double, Double) -> Exp Double -> Exp (Double, Double, Double)
div3 (unlift -> (a,b,c)) x = lift (a / x, b / x, c / x)

myfor :: Int -> (Int -> a -> a) -> a -> a
myfor n f x | n P.== 0    = x
            | P.otherwise =  myfor (n-1) f (f (n-1) x)

liftTupf :: (Elt a, Elt b) => (Exp a -> Exp a) -> (Exp b -> Exp b) -> Exp (a,b) -> Exp (a,b)
liftTupf f g (unlift -> (a, b)) = lift (f a, g b)

--Maybe backpermute could be used aswel, but we need source values
padder :: Elt a => Acc (Matrix a) -> Exp (Int, Int) -> Exp (Int, Int) -> Exp a ->  Acc (Matrix a)
padder array pad_width_x pad_width_y constant_val =
    let
        (x0, x1) = unlift pad_width_x :: (Exp Int, Exp Int)
        (y0, y1) = unlift pad_width_x :: (Exp Int, Exp Int)
        Z :. m :. n = (unlift . shape) array :: ( Z :. Exp Int :. Exp Int)
        def = fill (index2 (m + y0 + y1) (n + x0 + x1)) constant_val

        indexer (unlift -> Z :. y :. x :: Z :. Exp Int :. Exp Int)
            = index2 (y + y0) (x + x0)
        result = permute const def indexer array
    in result

c0 :: Exp Int
c0 = constant 0

addMessage :: String -> a -> a
addMessage mes action =
    let io = unsafePerformIO (P.putStrLn mes)
    in P.seq action (P.seq io action)

-- Check if the coordinates are out of bounds, and sets them to zero otherwise
fixoutofbounds :: Elt a => Exp Int -> Exp Int -> Exp a -> Exp (Int, Int, a) -> Exp (Int, Int, a)
fixoutofbounds maxx maxy defv
                old@(unlift -> (x', y', _) ::(Exp Int,Exp Int,Exp a)) =
    let outofbounds = x' < c0 || y' < c0 || x' >= maxx || y' >= maxy
        -- If out of bounds, set all values to zero (after the offset alters it)
        defx = c0 :: Exp Int
        defy = c0 :: Exp Int
        def = lift (defx, defy, defv)
    in if outofbounds then def else old