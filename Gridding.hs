{-# language FlexibleContexts    #-}
{-# language RebindableSyntax    #-}
{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}

module Gridding where

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A
import Data.Array.Accelerate.Math.DFT.Centre              as A hiding (shift1D, shift2D, shift3D)
import Data.Array.Accelerate.Math.FFT                     as A 

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..), maybe)

import Control.Lens as L (_1, _2)

import Data.Maybe (fromJust, fromMaybe)

data Kwargs = Kwargs { convolutionKernel :: Maybe (Array DIM4 (Complex Double))
                     , wSteps :: Maybe Int
                     , kernelCache :: Maybe (Int -> Array DIM4 (Complex Double))
                     -- TODO: Make good type for kernel_fn
                     , kernelFunction :: Maybe (Double -> Int -> Array DIM4 (Complex Double))
                     , patHorShift :: Maybe Int
                     , patVerShift :: Maybe Int
                     , patTransMat :: Maybe (Acc (Matrix Double))
                     }

noArgs :: Kwargs
noArgs = Kwargs Nothing Nothing Nothing Nothing Nothing Nothing Nothing
-------------------------------
-- Gridding Accelerate code
type ImagingFunction = Double                                        -- (theta) Field of view size
                     -> Int                                          -- (lam) Grid size
                     -> Acc (Vector (Double, Double, Double))        -- (uvw) all the uvw baselines (coordinates) (lenght : n * (n-1))
                     -> Acc (Vector (Int, Int, Double, Double))      -- (src) (Antenna 1, Antenna 2, The time (in MJD UTC), Frequency (Hz)) (lenght : n * (n-1))
                     -> Acc (Vector (Complex Double))                -- (vis) visibility  (length n * (n-1))
                     -> Kwargs                                       -- (kwargs) The additional options that you can give
                     -> Acc (Matrix (Complex Double))

-- # Simple imaging
simple_imaging :: ImagingFunction
simple_imaging theta lam uvw src vis kwargs =
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
conv_imaging theta lam uvw src vis kwargs =
    let kv = fromJust $ convolutionKernel kwargs
        lamf = fromIntegral lam
        n = P.round(theta * lamf)
        ne = constant n

        p = map (`div3` constant lamf) uvw
        czero = constant 0
        guv = fill (index2 ne ne) czero
    in convgrid kv guv p vis

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

convgrid :: Array DIM4 (Complex Double)                  -- The oversampled convolution kernel
         -> Acc (Matrix (Complex Double))                -- Destination grid N x N of complex numbers
         -> Acc (Vector (Double, Double, Double))        -- The uvw baselines, but between -.5 and .5
         -> Acc (Vector (Complex Double))                -- The visiblities
         -> Acc (Matrix (Complex Double))
convgrid gcf a p v =
    let Z :. qpx :. _ :. gh :. gw = (arrayShape gcf) :: Z :. Int :. Int :. Int :. Int
        Z :. height :. width = unlift (shape a) :: Z :. Exp Int :. Exp Int
        coords = frac_coords (lift (height, width)) (constant qpx) p
        gpugcf = use gcf
        halfgh = constant $ gh `div` 2
        halfgw = constant $ gw `div` 2

        getComplex :: Int -> Int -> Exp (Int, Int, Int, Int) -> Exp (Int, Int, Complex Double)
        getComplex i j (unlift -> (x, xf, y, yf)::(Exp Int,Exp Int,Exp Int,Exp Int)) =
            lift ( x
                 , y
                 , gpugcf ! lift (Z :. yf :. xf :. constant i :. constant j) )
        fullrepli :: Int -> Int -> Acc (Matrix (Complex Double)) -> Acc (Matrix (Complex Double))
        fullrepli i j origin = 
            let temp = map (getComplex i j) coords
                (x, y, source) = unzip3 temp
                indexer id = lift (Z :. y ! id - halfgh + constant i :. x ! id - halfgw + constant j)
                newV = zipWith (*) source v
            in permute (+) origin indexer newV
        
    in myfor gh (\i -> myfor gw (fullrepli i)) a


-- Imaging with a w kernel
w_cache_imaging :: ImagingFunction
w_cache_imaging theta lam uvw src vis kwargs@Kwargs{wSteps = wstep, kernelCache = kernel_cache, kernelFunction = kernel_fn} =
    let lamf = fromIntegral lam
        n = P.round(theta * lamf)
        ne = constant n

        p = map (`div3` constant lamf) uvw
        czero = constant 0
        guv = fill (index2 ne ne) czero :: Acc (Matrix (Complex Double))
    in undefined
----------------------------------
-- Processing the imaging functions
do_imaging :: Double                                 -- Field of view size
           -> Int                                    -- Grid size
           -> Vector (Double, Double, Double)        -- all the uvw baselines (coordinates) (lenght : n * (n-1))
           -> Vector (Int, Int, Double, Double)      -- (Antenna 1, Antenna 2, The time (in MJD UTC), Frequency (Hz)) (lenght : n * (n-1))
           -> Vector (Complex Double)                -- visibility  (length n * (n-1))
           -> ImagingFunction
           -> Kwargs
           -> Acc (Matrix Double,Matrix Double, Scalar Double)
do_imaging theta lam uvw src vis imgfn kwargs = 
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
        cdrt = imgfn theta lam uvw1 src0 (zipWith (*) wt vis1) kwargs
        drt = (map real . ifft . make_grid_hermitian) cdrt
        -- Make point spread function
        c = imgfn theta lam uvw1 src0 wt kwargs
        psf = (map real . ifft . make_grid_hermitian) c
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
w_kernel :: Exp Double                           -- Field of view
         -> Exp Double                           -- Baseline distance to the projection plane
         -> Exp Int                              -- Far field size
         -> Exp Int                              -- Size of convolution
         -> Exp Int                              -- Oversampling
         -> Kwargs                               -- Additional kernel arguments               
         -> Acc (Array DIM4 (Complex Double))    -- [Qpx,Qpx,s,s] shaped oversampled convolution kernels             
w_kernel theta w npixFF npixKern qpx kwargs =
    let
        (l, m) = unlift $ kernel_coordinates npixFF theta kwargs :: (Acc (Matrix Double), Acc (Matrix Double))
        kern = w_kernel_function l m w
    in kernel_oversample kern npixFF qpx npixKern

kernel_coordinates :: Exp Int -> Exp Double -> Kwargs -> Acc ((Matrix Double), (Matrix Double))
kernel_coordinates n theta kwargs =
    let
        dl = (constant . fromIntegral . fromMaybe 0 . patHorShift) kwargs
        dm = (constant . fromIntegral . fromMaybe 0 . patVerShift) kwargs
        (l,m) = unlift $ coordinates2 (n) :: (Acc (Matrix Double), Acc (Matrix Double))-- * theta 
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
                  -> Exp Double                       -- Baseline distance to the projection plane
                  -> Acc (Matrix (Complex Double))    --N x N array with the far field
w_kernel_function l m w = 
    let r2 = zipWith (\x y -> x**x + y*y) l m
    {-  dl and dm seem to be always zero anyway in the reference code
        dl = 0
        dm = 0
        ph1 = zipWith3 (\r2' l' m' -> 1 - sqrt (1 - r2') - dl * l' - dm * m') r2 l m
        ph2 = map (\r2' -> 1 - sqrt (1 - r2')) r2
        ph  = if dl == 0 && dm == 0 then ph2 else ph1
    -}
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
        pad_width = lift ((n `div` 2)-(n0 `div` 2), (n+1 `div` 2)-(n0+1 `div` 2))
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
        qpx2 = fromIntegral (qpx * qpx)
    in map (*qpx2) w_kern

----------------------
-- Fourier transformations
fft :: Acc (Matrix (Complex Double)) -> Acc (Matrix (Complex Double))
fft = shift2D . fft2D Forward . ishift2D

ifft :: Acc (Matrix (Complex Double)) -> Acc (Matrix (Complex Double))
ifft = shift2D . fft2D Inverse . ishift2D

shift1D :: Elt e => Acc (Vector e) -> Acc (Vector e)
shift1D arr = backpermute sh p arr  
      where
        sh      = shape arr
        n       = indexHead sh
        --
        shift   = (n `quot` 2) + boolToInt (odd n)
        roll i  = (i+shift) `rem` n
        p       = ilift1 roll

ishift1D :: Elt e => Acc (Vector e) -> Acc (Vector e)
ishift1D arr = backpermute sh p arr  
      where
        sh      = shape arr
        n       = indexHead sh
        --
        shift   = (n `quot` 2)-- + boolToInt (odd n)
        roll i  = (i+shift) `rem` n
        p       = ilift1 roll

shift2D :: Elt e => Acc (Array DIM2 e) -> Acc (Array DIM2 e)
shift2D arr
  = backpermute sh p arr
  where
    sh      = shape arr
    Z :. h :. w = unlift sh
    --
    shifth = (h `quot` 2) + boolToInt (odd h)
    shiftw = (w `quot` 2) + boolToInt (odd w)

    p ix
      = let Z:.y:.x = unlift ix :: Z :. Exp Int :. Exp Int
        in index2 ((y + shifth) `rem` h)
                  ((x + shiftw) `rem` w)

ishift2D :: Elt e => Acc (Array DIM2 e) -> Acc (Array DIM2 e)
ishift2D arr
  = backpermute sh p arr
  where
    sh      = shape arr
    Z :. h :. w = unlift sh
    --
    shifth = (h `quot` 2)
    shiftw = (w `quot` 2)

    p ix
      = let Z:.y:.x = unlift ix :: Z :. Exp Int :. Exp Int
        in index2 ((y + shifth) `rem` h)
                  ((x + shiftw) `rem` w)

shift3D :: Elt e => Acc (Array DIM3 e) -> Acc (Array DIM3 e)
shift3D arr
  = backpermute sh p arr
  where
    sh      = shape arr
    Z :. d :. h :. w = unlift sh
    --
    shiftd = (d `quot` 2) + boolToInt (odd d)
    shifth = (h `quot` 2) + boolToInt (odd h)
    shiftw = (w `quot` 2) + boolToInt (odd w)

    p ix
      = let Z:.z:.y:.x = unlift ix :: Z :. Exp Int :. Exp Int :. Exp Int
        in index3 ((z + shiftd) `rem` d)
                  ((y + shifth) `rem` h)
                  ((x + shiftw) `rem` w)

ishift3D :: Elt e => Acc (Array DIM3 e) -> Acc (Array DIM3 e)
ishift3D arr
  = backpermute sh p arr
  where
    sh      = shape arr
    Z :. d :. h :. w = unlift sh
    --
    shiftd = (d `quot` 2)
    shifth = (h `quot` 2)
    shiftw = (w `quot` 2)

    p ix
      = let Z:.z:.y:.x = unlift ix :: Z :. Exp Int :. Exp Int :. Exp Int
        in index3 ((z + shiftd) `rem` d)
                  ((y + shifth) `rem` h)
                  ((x + shiftw) `rem` w)

------------------------
-- Helper functions
div3 :: Exp (Double, Double, Double) -> Exp Double -> Exp (Double, Double, Double)
div3 (unlift -> (a,b,c)) x = lift (a / x, b / x, c / x)

myfor :: Int -> (Int -> a -> a) -> a -> a
myfor n f x | n P.== 0  = x
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