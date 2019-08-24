{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE RebindableSyntax    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators       #-}
{-# LANGUAGE ViewPatterns        #-}

module Gridding where

import           FFT
import           Types
import           Utility

import           Data.Array.Accelerate                 as A hiding (fromInteger,
                                                             fromIntegral,
                                                             fromRational)
import qualified Data.Array.Accelerate                 as A (fromInteger,
                                                             fromIntegral,
                                                             fromRational)
import           Data.Array.Accelerate.Data.Complex    as A
import           Data.Array.Accelerate.Math.DFT.Centre as A
import           Data.Array.Accelerate.Math.FFT        as A

import           Data.Array.Accelerate.Interpreter     as I
import           Data.Array.Accelerate.LLVM.Native     as CPU

import           Prelude                               as P (IO, Maybe (..),
                                                             String,
                                                             fromInteger,
                                                             fromIntegral,
                                                             fromRational,
                                                             maybe, return,
                                                             (>>), (>>=))
import qualified Prelude                               as P

import           Control.Lens                          as L (_1, _2, _3, _4, _5)
import           Data.Maybe                            (fromJust, fromMaybe)

import           Debug.Trace

data KernelOptions = KernelOptions
    { patHorShift :: Maybe Int
    , patVerShift :: Maybe Int
    , patTransMat :: Maybe (Acc (Matrix F))
    , wstep       :: Maybe Int
    , qpx         :: Maybe Int              -- The oversampling of the convolution kernel
    , npixFF      :: Maybe Int
    , npixKern    :: Maybe Int         -- Kernel size (if it is same size in x and y)
    }

data OtherImagingArgs = OtherImagingArgs
                      { convolutionKernel :: Maybe Kernel
                      , akernels :: Maybe (Array DIM3 Visibility)
                      , wkernels :: Maybe ((Array DIM5 Visibility, Vector BaseLine))
                      , kernelCache :: Maybe KernelF
                      , kernelFunction :: Maybe KernelF
                      }

type KernelF = Exp F                                 -- Field of view
             -> Acc (Scalar BaseLine)                -- Baseline distance to the projection plane
             -> Maybe (Exp Antenna)                  -- The id of the antenna 1
             -> Maybe (Exp Antenna)                  -- The id of the antenna 2
             -> Maybe (Exp Time)                     -- The time
             -> Maybe (Exp Frequency)                -- The frequency
             -> KernelOptions                        -- Additional kernel arguments
             -> Acc Kernel                           -- [Qpx,Qpx,s,s] shaped oversampled convolution kernels

type WKernelF = Exp F                                 -- Field of view
              -> Acc (Scalar BaseLine)                -- Baseline distance to the projection plane
              -> KernelOptions                        -- Additional kernel arguments
              -> Acc Kernel                           -- [Qpx,Qpx,s,s] shaped oversampled convolution kernels

type AKernelF = Exp F                                 -- Field of view
              -> Exp Antenna                          -- The id of the antenna
              -> Exp Time                             -- The time
              -> Exp Frequency                        -- The frequency
              -> Acc (Matrix Visibility)              -- The a-kernel


noArgs :: KernelOptions
noArgs = KernelOptions Nothing Nothing Nothing Nothing Nothing Nothing Nothing

noOtherArgs :: OtherImagingArgs
noOtherArgs = OtherImagingArgs Nothing Nothing Nothing Nothing Nothing
-------------------------------
-- Gridding Accelerate code
type ImagingFunction = F                                                 -- (theta) Field of view size
                     -> Int                                              -- (lam) Grid size
                     -> Acc (Vector BaseLines)                           -- (uvw) all the uvw baselines (coordinates) (lenght : n * (n-1))
                     -> Acc (Vector (Antenna, Antenna, Time, Frequency)) -- (src) (Antenna 1, Antenna 2, The time (in MJD UTC), Frequency (Hz)) (lenght : n * (n-1))
                     -> Acc (Vector Visibility)                          -- (vis) visibility  (length n * (n-1))
                     -> Acc (Matrix Visibility)

-- # Simple imaging
simple_imaging :: ImagingFunction
simple_imaging theta lam uvw src vis =
    let lamf = fromIntegral lam
        n = P.round(theta * lamf)
        ne = constant n

        p = map (`div3` constant lamf) uvw
        czero = constant 0
        guv = fill (index2 ne ne) czero
    in grid guv p vis

grid :: Acc (Matrix Visibility)                    -- Destination grid N x N of complex numbers
     -> Acc (Vector BaseLines)                     -- The uvw baselines, but between -.5 and .5
     -> Acc (Vector Visibility)                    -- The visiblities
     -> Acc (Matrix Visibility)
grid a p v = permute (+) a (\i -> xy ! i) v
    where
        Z :. n :. _ = unlift (shape a) :: Z :. Exp Int :. Exp Int
        halfn = n `div` 2
        nf = A.fromIntegral n :: Exp F
        xy = map gridxy p :: Acc (Vector DIM2)

        -- Note the order, python specified as a[y, x] += v
        -- We do the same. This is correct, first y then x in the index2 function
        gridxy :: Exp BaseLines -> Exp DIM2
        gridxy (unlift -> (x, y, _ :: Exp F)) = index2 (toGridCell y) (toGridCell x)

        toGridCell :: Exp F -> Exp Int
        toGridCell f = halfn + floor (0.5 + nf * f)

-- # Normal imaging
conv_imaging :: Kernel -> ImagingFunction
conv_imaging kv theta lam uvw src vis =
    let lamf = fromIntegral lam
        n = P.round(theta * lamf)
        ne = constant n

        p = map (`div3` constant lamf) uvw
        czero = constant 0
        guv = fill (index2 ne ne) czero
    in convgrid (use kv) guv p vis

frac_coord :: Exp Int                                   -- The height/width
           -> Exp Int                                   -- Oversampling factor
           -> Acc (Vector BaseLine)                     -- a uvw baseline, but between -.5 and .5
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
            -> Acc (Vector BaseLines)                   -- The uvw baselines, but between -.5 and .5
            -> Acc (Vector (Int, Int, Int, Int))
-- NOTE we should give it in height, width order.
frac_coords (unlift -> (h, w) ) qpx p =
    let (u, v, _) = unzip3 p
        (x, xf)   = unzip (frac_coord w qpx u)
        (y, yf)   = unzip (frac_coord h qpx v)
    in zip4 x xf y yf

convgrid :: Acc (Array DIM4 Visibility)        -- The oversampled convolution kernel
         -> Acc (Matrix Visibility)            -- Destination grid N x N of complex numbers
         -> Acc (Vector BaseLines)             -- The uvw baselines, but between -.5 and .5
         -> Acc (Vector Visibility)            -- The visiblities
         -> Acc (Matrix Visibility)
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

        getComplexAndAddOffset :: Exp DIM3 -> Exp (Int, Int, Int, Int, Visibility) -> Exp (Int, Int, Visibility)
        getComplexAndAddOffset
            (unlift . unindex3 -> (_, i, j) ::(Exp Int,Exp Int,Exp Int))
            (unlift -> (x, xf, y, yf, vis)::(Exp Int,Exp Int,Exp Int,Exp Int,Exp Visibility)) =
            lift ( x + j
                 , y + i
                 , vis * gcf ! lift (Z :. yf :. xf :. i :. j) )
    in permute (+) a indexer val

convgrid2 :: Acc (Array DIM5 Visibility)        -- The oversampled convolution kernel
          -> Acc (Matrix Visibility)            -- Destination grid N x N of complex numbers
          -> Acc (Vector BaseLines)             -- The uvw baselines, but between -.5 and .5
          -> Acc (Vector Int)                   -- *DIF* The wbin index of the convolution kernel
          -> Acc (Vector Visibility)            -- The visiblities
          -> Acc (Matrix Visibility)
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

        getComplexAndAddOffset :: Exp DIM3 -> Exp (Int, Int, Int, Int, Int, Visibility) -> Exp (Int, Int, Visibility)
        getComplexAndAddOffset
            (unlift . unindex3 -> (_, i, j) ::(Exp Int,Exp Int,Exp Int))
            (unlift -> (wbin, x, xf, y, yf, vis)::(Exp Int,Exp Int,Exp Int,Exp Int,Exp Int,Exp Visibility)) =
            lift ( x + j
                 , y + i
                 , vis * gcf ! lift (Z :. wbin :. yf :. xf :. i :. j) )
    in permute (+) a indexer val

convgrid3 ::
             Acc (Array DIM5 Visibility)        -- The oversampled convolution w-kernel
          -> Acc (Array DIM3 Visibility)        -- The a-kernels
          -> Acc (Matrix Visibility)            -- Destination grid N x N of complex numbers
          -> Acc (Vector BaseLines)             -- The uvw baselines, but between -.5 and .5
          -> Acc (Vector (Int, Int, Int))       -- *DIF* The wbin index of the convolution kernel and the index of the a-kernels
          -> Acc (Vector Visibility)            -- The visiblities
          -> Acc (Matrix Visibility)
convgrid3 wkerns akerns a p index v =
    let Z :. _ :. qpx :. _ :. gh :. gw = unlift (shape wkerns) :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int :. Exp Int
        Z :. height :. width = unlift (shape a) :: Z :. Exp Int :. Exp Int
        Z :. n = unlift (shape v) :: Z :. Exp Int
        halfgh = gh `div` 2
        halfgw = gw `div` 2

        --Gather the coords, make them into ints and fractions and zip with visibility
        coords = frac_coords (lift (height, width)) qpx p
        (cx, cxf, cy, cyf) = unzip4 coords
        dat = zip5 cx cxf cy cyf v

        processer :: Exp Int -> Acc (Matrix Visibility) -> Acc (Matrix Visibility)
        processer n aa =
            let id = unit $ index !! n
                dat = unit $ dat !! n
            in processOne id dat wkerns akerns aa
    in afor n processer a

processOne :: Acc (Scalar (Int, Int, Int)) -> Acc (Scalar(Int, Int, Int, Int, Visibility)) -> Acc (Array DIM5 Visibility) -> Acc (Array DIM3 Visibility) -> Acc (Matrix Visibility) -> Acc (Matrix Visibility)
processOne
    (unlift . the -> (wbin, a1index, a2index) :: (Exp Int, Exp Int, Exp Int))
    (unlift . the -> (x, xf, y, yf, vis)::(Exp Int,Exp Int,Exp Int,Exp Int,Exp Visibility))
    wkerns akerns a =
        let Z :. _ :. _ :. _ :. gh :. gw = unlift (shape wkerns) :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int :. Exp Int
            Z :. height :. width = unlift (shape a) :: Z :. Exp Int :. Exp Int
            halfgh = gh `div` 2
            halfgw = gw `div` 2
            cc0 = (lift (constant 0.0 :+ constant 0.0))

            --Get the right aw kernel
            a1 = slice akerns (lift (Z :. a1index :. All :. All))
            a2 = slice akerns (lift (Z :. a2index :. All :. All))
            w  = slice wkerns (lift (Z :. wbin :. All :. All :. All :. All))
            -- NOTE, the conjugate normally happens in imaging function, but it is convenient to do it here.
            awkern = map conjugate $ aw_kernel_fn2 yf xf w a1 a2

            --Start at the right coordinates
            startx = x - halfgw
            starty = y - halfgh
            vals = fill (shape awkern) (lift (startx, starty, vis))
            mappedvals = imap getComplexAndAddOffset vals
            --Fix the bounds
            fixer = fixoutofbounds width height cc0
            fixedBounds = map fixer mappedvals
            (xs,ys, val) = unzip3 fixedBounds

            indexer id =
                let y' = ys ! id
                    x' = xs ! id
                in index2 y' x'

            getComplexAndAddOffset :: Exp DIM2 -> Exp (Int, Int, Visibility) -> Exp (Int, Int, Visibility)
            getComplexAndAddOffset (unlift . unindex2 -> (i, j)::(Exp Int,Exp Int))
                (unlift -> (x, y, vis)::(Exp Int,Exp Int,Exp Visibility)) =
                    lift ( x + j
                         , y + i
                         , vis * awkern ! lift (Z :. i :. j) )
        in permute (+) a indexer val

convgridSeq' :: (Mode -> Acc (Matrix Visibility) -> Acc (Matrix  Visibility))
             -> Acc (Array DIM5 Visibility)        -- The oversampled convolution w-kernel
             -> Acc (Array DIM3 Visibility)        -- The a-kernels
             -> Acc (Matrix Visibility)            -- Destination grid N x N of complex numbers
             -> Acc (Vector BaseLines)             -- The uvw baselines, but between -.5 and .5
             -> Acc (Vector (Int, Int, Int))       -- *DIF* The wbin index of the convolution kernel and the index of the a-kernels
             -> Acc (Vector Visibility)            -- The visiblities
             -> Acc (Matrix Visibility)
convgridSeq' myfft wkerns akerns a p index v =
    let Z :. _ :. qpx :. _ :. gh :. gw = unlift (shape wkerns) :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int :. Exp Int
        Z :. height :. width = unlift (shape a) :: Z :. Exp Int :. Exp Int
        halfgh = gh `div` 2
        halfgw = gw `div` 2

        --Gather the coords, make them into ints and fractions and zip with visibility
        coords             = frac_coords (lift (height, width)) qpx p
        (cx, cxf, cy, cyf) = unzip4 coords
        (wbin, a1, a2)     = A.unzip3 index
        visIndex           = A.zip6 wbin a1 a2 cxf cyf v

        visSeq     = toSeqInner visIndex
        coordsSeq  = toSeqInner $ A.zip cx cy
        visKernSeq = mapSeq (processOneSeq' myfft wkerns akerns) visSeq

        addCoords :: Acc (Matrix Visibility) -> Acc (Scalar (Int, Int)) -> Acc (Matrix (Int, Int, Visibility))
        addCoords vis xy = let
            (x,y) = unlift . the $ xy :: (Exp Int, Exp Int)
            indexmapper (unlift . unindex2 -> (i, j)::(Exp Int,Exp Int)) vis =
                    lift ( x + i - halfgw, y + j - halfgh, vis)
            in imap indexmapper vis

        visAndCoordSeq = zipWithSeq addCoords visKernSeq coordsSeq
        res = collect . elements $ visAndCoordSeq
        (xs, ys, resvis) = A.unzip3 res

        indexer id = let y' = ys ! id
                         x' = xs ! id
                     in index2 y' x'
    in permute (+) a indexer resvis

convgridSeq :: Acc (Array DIM5 Visibility)        -- The oversampled convolution w-kernel
            -> Acc (Array DIM3 Visibility)        -- The a-kernels
            -> Acc (Matrix Visibility)            -- Destination grid N x N of complex numbers
            -> Acc (Vector BaseLines)             -- The uvw baselines, but between -.5 and .5
            -> Acc (Vector (Int, Int, Int))       -- *DIF* The wbin index of the convolution kernel and the index of the a-kernels
            -> Acc (Vector Visibility)            -- The visiblities
            -> Acc (Matrix Visibility)
convgridSeq = convgridSeq' myfft2D

processOneSeq' :: (Mode -> Acc (Matrix Visibility) -> Acc (Matrix Visibility)) 
               -> Acc (Array DIM5 Visibility) -> Acc (Array DIM3 Visibility) -> Acc (Scalar (Int, Int, Int, Int, Int, Visibility)) ->  Acc (Matrix Visibility)
processOneSeq' myfft wkerns akerns
    (unlift . the -> (wbin, a1index, a2index, xf, yf, vis) :: (Exp Int, Exp Int, Exp Int, Exp Int, Exp Int, Exp Visibility)) =
        let Z :. _ :. _ :. _ :. gh :. gw = unlift (shape wkerns) :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int :. Exp Int
            halfgh = gh `div` 2
            halfgw = gw `div` 2

            --Get the right a and w kernels
            a1 = slice akerns (lift (Z :. a1index :. All :. All))
            a2 = slice akerns (lift (Z :. a2index :. All :. All))
            w  = slice wkerns (lift (Z :. wbin :. yf :. xf :. All :. All))
            --Convolve them
            -- NOTE, the conjugate normally happens in imaging function, but it is convenient to do it here.
            akern = convolve2d' myfft a1 a2
            awkern = map conjugate $ convolve2d' myfft w akern
            --awkern = convolve3_2d a1 a2 w

            -- Let the visibility have the same dimensions as the aw-kernel
            allvis = replicate ( lift ( Z :. gh :. gw)) (unit vis)
        in zipWith (*) allvis awkern

processOneSeq ::  Acc (Array DIM5 Visibility) -> Acc (Array DIM3 Visibility) -> Acc (Scalar (Int, Int, Int, Int, Int, Visibility)) ->  Acc (Matrix Visibility)
processOneSeq = processOneSeq' myfft2D
{-
-- Imaging with a w kernel TODO: It differs from the normal implementation atm, it will only work with w-kernels, not aw-kernels
w_cache_imaging :: KernelOptions -> OtherImagingArgs -> ImagingFunction
w_cache_imaging kernops@KernelOptions{wstep = wstep'}
  otargs@OtherImagingArgs{kernelCache = kernel_cache', kernelFunction = kernel_fn'}
  theta lam uvw src vis
  =
    let
        -- They use a LRU cache for the function, to memoize the results. Maybe usefull. Not sure yet.
        -- So maybe use LRU cache implementation (lrucache package) or memoization (memoize package)
        -- For now we just call the function.
        default_kernelf theta w _ _ _ _ kernops = w_kernel theta w kernops
        kernel_fn    = fromMaybe default_kernelf kernel_fn'
        kernel_cache = fromMaybe kernel_fn kernel_cache'
        only_wcache theta w kernops = kernel_cache theta w Nothing Nothing Nothing Nothing kernops
        wstep_ = fromMaybe 2000 wstep'
        wstep = constant $ wstep_

        lamf = fromIntegral lam
        n = P.round(theta * lamf)
        ne = constant n

        p = map (`div3` constant lamf) uvw
        czero = constant 0
        guv = fill (index2 ne ne) czero :: Acc (Matrix Visibility)


        -- We just make them all here, that's easiest for now
        (_,_,w) = unzip3 uvw
        roundedw = map (\w' -> wstep * round ( w' / A.fromIntegral wstep)) w
        maxw = the $ maximum roundedw
        minw = the $ minimum roundedw
        steps = ((maxw - minw) `div` wstep ) + 1
        (cpumin, cpumax, cpusteps) = CPU.run (unit $ lift (minw, maxw, steps)) `indexArray` Z

        wbins = map (\rw' -> (rw' - constant cpumin) `div` wstep) roundedw

        makeWKernel :: Int -> Acc (Array DIM5 Visibility)
        makeWKernel i = let wbin = fromList Z $ (fromIntegral $ i * wstep_ + cpumin) : []
                        in use $ compiledMakeWKernel wbin

        compiledMakeWKernel = CPU.runN makeWKernel'

        makeWKernel' :: Acc (Scalar F) -> Acc (Array DIM5 Visibility)
        makeWKernel' wbin = make5D . map conjugate $ only_wcache (constant theta) wbin kernops


        make5D mat = let (Z :. yf :. xf :. y :. x) =  (unlift . shape) mat :: (Z :.Exp Int :. Exp Int :. Exp Int :. Exp Int)
                         newsh = lift (Z :. constant 1 :. yf :. xf :. y :. x) :: Exp DIM5
                     in reshape newsh mat

        thekernels = CPU.runN $ compute $ myfor (cpusteps - 1) (\i old -> concatOn _5 (makeWKernel i) old) (makeWKernel (cpusteps - 1))
    in convgrid2 (use thekernels) guv p wbins vis
-}

-- Imaging with aw-caches
aw_imaging' :: (Mode -> Acc (Matrix Visibility) -> Acc (Matrix Visibility))
            -> KernelOptions -> OtherImagingArgs -> F -> Int
            -> Acc WKernels
            -> Acc (Vector F)
            -> Acc AKernels
            -> Acc (Vector BaseLines)                           -- (uvw) all the uvw baselines (coordinates) (lenght : n * (n-1))
            -> Acc (Vector (Antenna, Antenna, Time, Frequency)) -- (src) (Antenna 1, Antenna 2, The time (in MJD UTC), Frequency (Hz)) (lenght : n * (n-1))
            -> Acc (Vector Visibility)                          -- (vis) visibility  (length n * (n-1))
            -> Acc (Matrix Visibility)
aw_imaging' myfft kernops@KernelOptions{wstep = wstep', qpx = qpx'}
  otargs@OtherImagingArgs{akernels = akernels'
                         ,wkernels = wkernels'} theta lam
        wkernels wbins akernels uvw src vis =
    let
        lamf = fromIntegral lam
        n = constant . P.round $ theta * lamf

        p = map (`div3` constant lamf) uvw
        czero = constant 0
        guv = fill (index2 n n) czero :: Acc (Matrix Visibility)

        -- We just make them all here, that's easiest for now
        (_,_,w) = unzip3 uvw
        closestw = map (findClosest wbins) w
        (a1, a2, _, _) = unzip4 src
        index = zipWith3 (\x y z -> lift (x, A.fromIntegral y, A.fromIntegral z)) closestw a1 a2
        --Normally we conjugate the kernel here, but we do it later on in the convgrid3 (or processOne) function somewhere
    in convgridSeq' myfft wkernels akernels guv p index vis

aw_imaging :: KernelOptions -> OtherImagingArgs -> F -> Int
           -> Acc WKernels
           -> Acc (Vector F)
           -> Acc AKernels
           -> Acc (Vector BaseLines)                           -- (uvw) all the uvw baselines (coordinates) (lenght : n * (n-1))
           -> Acc (Vector (Antenna, Antenna, Time, Frequency)) -- (src) (Antenna 1, Antenna 2, The time (in MJD UTC), Frequency (Hz)) (lenght : n * (n-1))
           -> Acc (Vector Visibility)                          -- (vis) visibility  (length n * (n-1))
           -> Acc (Matrix Visibility)
aw_imaging = aw_imaging' myfft2D
----------------------------------
-- Processing the imaging functions
do_imaging :: F                                -- Field of view size
           -> Int                              -- Grid size
           -> Matrix BaseLine                  -- all the uvw baselines (coordinates) (dims : (n * (n-1), 3))
           -> Vector Antenna                   -- Antenna 1 (lenght : n * (n-1))
           -> Vector Antenna                   -- Antenna 2
           -> Vector Time                      -- The time (in MJD UTC)
           -> Frequency                        -- Frequency (Hz)
           -> Vector Visibility                -- visibility  (length n * (n-1))
           -> ImagingFunction
           -> Acc (Image,Image, Scalar F)
do_imaging theta lam uvw a1 a2 t f vis imgfn =
    let
        -- TODO: if src == None: src = numpy.ndarray((len(vis), 0))
        len = arrayShape vis
        uvwmat = use uvw
        u = slice uvwmat (constant (Z :. All :. (0 :: Int)))
        v = slice uvwmat (constant (Z :. All :. (1 :: Int)))
        w = slice uvwmat (constant (Z :. All :. (2 :: Int)))
        uvw0 = zip3 u v w
        vis0 = use vis
        src0 = zip4 (use a1) (use a2) (use t) (fill (lift len) (constant f))
        -- Mirror baselines such that v>0
        (uvw1,vis1) = unzip $ mirror_uvw uvw0 vis0

        -- Determine weights
        ones = fill (lift len) 1
        wt = doweight theta lam uvw1 ones

        -- Make image
        cdrt = imgfn theta lam uvw1 src0 (zipWith (*) wt vis1)
        drt = (map real . ifft . make_grid_hermitian) cdrt
        -- Make point spread function
        c = imgfn theta lam uvw1 src0 wt
        psf =  (map real . ifft . make_grid_hermitian) c
        -- Normalise
        pmaxv = (maximum . flatten) psf
        pmax = the pmaxv
        res0 = map (/ pmax) drt
        res1 = map (/ pmax) psf
        res2 = pmaxv
    in lift (res0, res1, res2)

mirror_uvw :: Acc (Vector BaseLines)                              -- all the uvw baselines (coordinates) (lenght : n * (n-1))
           -> Acc (Vector Visibility)                             -- visibility  (length n * (n-1))
           -> Acc (Vector (BaseLines, Visibility))
mirror_uvw uvw vis =
    let
        (u,v,w) = unzip3 uvw
        uvwvis = zip4 u v w vis
        mirrorv o@(unlift -> (u,v,w,vis) :: (Exp BaseLine, Exp BaseLine, Exp BaseLine, Exp Visibility)) =
            if v < 0
                then lift ((-u, -v, -w), conjugate vis)
                else lift ((u,v,w),vis)
    in map mirrorv uvwvis

doweight :: F                                                  -- Field of view size
         -> Int                                                -- Grid size
         -> Acc (Vector BaseLines)                             -- all the uvw baselines (coordinates) (lenght : n * (n-1))
         -> Acc (Vector Visibility)                            -- Array of ones (same size as visibility  (length n * (n-1)))
         -> Acc (Vector Visibility)
doweight theta lam p v =
    let
        n = P.round(theta * lamf)
        lamf = fromIntegral lam
        ne = constant n :: Exp Int
        lamfe = constant lamf
        gw = fill (index2 ne ne) 0 :: Acc (Matrix F)
        coords = frac_coords (lift (ne, ne)) 1 (map (`div3` lamfe) p)
        (x,_,y,_) = unzip4 coords
        ones  = fill (shape x) 1 :: Acc (Vector F)
        xyindex = zipWith index2 y x
        weights = permute (+) gw (\id -> xyindex ! id) ones

        newv = imap (\id v -> v / lift (weights ! (xyindex ! id) :+ 0)) v
    in newv

make_grid_hermitian :: Acc (Matrix Visibility) -> Acc (Matrix Visibility)
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
        oddReshape = (reverse2 . reverse1) guv
        evenReshape = backpermute sh indexer guv

        oddConj = map conjugate oddReshape
        evenConj = imap mapper evenReshape

        reshapedConj = if even n then evenConj else oddConj
    in zipWith (+) guv reshapedConj


----------------------
-- Kernels
w_kernel :: WKernelF
w_kernel theta w
    kernops@KernelOptions{npixFF = npixFF', npixKern = npixKern', qpx = qpx'} =
    let
        npixFF = constant $ fromJust npixFF'
        npixKern = constant $ fromJust npixKern'
        qpx = constant $ fromJust qpx'
        (l, m) = unlift $ kernel_coordinates npixFF theta kernops :: (Acc (Matrix BaseLine), Acc (Matrix BaseLine))
        kern = w_kernel_function l m w
    in kernel_oversample kern npixFF qpx npixKern

kernel_coordinates :: Exp Int -> Exp F -> KernelOptions -> Acc ((Matrix BaseLine), (Matrix BaseLine))
kernel_coordinates n theta kernops =
    let
        dl = (constant . fromIntegral . fromMaybe 0 . patHorShift) kernops
        dm = (constant . fromIntegral . fromMaybe 0 . patVerShift) kernops
        (l,m) = unlift $ coordinates2 n :: (Acc (Matrix BaseLine), Acc (Matrix BaseLine))
        lm1 = map (liftTupf (*theta) (*theta)) (zip l m)
        lm2 = case patTransMat kernops of
            Nothing -> lm1
            Just t  -> map (\pair -> let (x, y) = unlift pair in lift
                ( (t ! index2 0 0) * x + (t ! index2 1 0) * y
                , (t ! index2 0 1) * x + (t ! index2 1 1) * y ) ) lm1

        lm3 = map (liftTupf (+dl) (+dm)) lm2
    in lift $ unzip lm3

coordinates2 :: Exp Int -> Acc (Matrix BaseLine, Matrix BaseLine)
coordinates2 n =
    let
        n2 = n `div` 2
        sh = index1 n
        sh2 = index2 n n
        start = A.fromIntegral (-n2) * step
        step = 1 / A.fromIntegral n
        base = enumFromStepN sh start step :: Acc (Vector BaseLine)
        samecolumns = replicate (lift (Z :. n :. All)) base
        samerows = replicate (lift (Z :. All :. n)) base
    in lift (samecolumns, samerows)


w_kernel_function :: Acc (Matrix BaseLine)              -- Horizontal image coordinates
                  -> Acc (Matrix BaseLine)              -- Vertical image coordinates
                  -> Acc (Scalar BaseLine)              -- Baseline distance to the projection plane
                  -> Acc (Matrix Visibility)    -- N x N array with the far field
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

kernel_oversample :: Acc (Matrix Visibility)             -- Far field pattern
                  -> Exp Int                             -- Image size without oversampling
                  -> Exp Int                             -- Factor to oversample by -- there will be Qpx x Qpx convolution arl
                  -> Exp Int                             -- Size of convolution function to extract
                  -> Acc Kernel                          -- array of shape [ov, ou, v, u], e.g. with sub-pixel
kernel_oversample ff n qpx s =
    let
        -- Pad the far field to the required pixel size
        padff = pad_mid ff (n*qpx)
        -- Obtain oversampled uv-grid
        af = ifft padff
    in extract_oversampled af qpx s

aw_kernel_fn :: Exp Int
             -> Exp Int
             -> AKernelF
             -> Maybe WKernelF
             -> KernelF
aw_kernel_fn yf xf a_kernel_fn w_kernel_fn' theta w a1' a2' t' f'
    kernops@KernelOptions{qpx = qpx', npixKern = n'} =
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

        wkern = w_kernel_fn theta w kernops

        getkern :: Elt a => Exp Int -> Exp Int -> Acc (Array DIM4 a) -> Acc (Matrix a)
        getkern yf xf a = slice a (lift (Z :. yf :. xf :. All :. All))

        newd = constant $ Z :. 1 :. 1 :. n :. n

        awkernf =  reshape newd $ convolve2d akern (getkern yf xf wkern)
    in awkernf

aw_kernel_fn2 :: Exp Int
              -> Exp Int
              -> Acc WKernel
              -> Acc AKernel
              -> Acc AKernel
              -> Acc AWKernel
aw_kernel_fn2 yf xf wkern a1kern a2kern =
    let
        akern = convolve2d a1kern a2kern

        getkern :: Elt a => Exp Int -> Exp Int -> Acc (Array DIM4 a) -> Acc (Matrix a)
        getkern yf xf a = slice a (lift (Z :. yf :. xf :. All :. All))

        awkern =  convolve2d akern (getkern yf xf wkern)
    in awkern

-- the kernels for Aw-gridding will normally be chosen to /not/ overflow the borders
-- Thus we chose the simpler method
convolve2dO :: Acc (Matrix Visibility) -> Acc (Matrix Visibility) -> Acc (Matrix Visibility)
convolve2dO a1 a2 =
    let
        (Z :. n_ :. _) = (unlift . shape) a1 :: Z :. Exp Int :. Exp Int
        pw = A.ceiling (logBase 2 (A.fromIntegral n_) :: Exp F) :: Exp Int
        n = 2 ^ pw :: Exp Int
        n2 = A.fromIntegral $ n * n

        a1fft = myfft2D Inverse . ishift2D $ pad_mid a1 n
        a2fft = myfft2D Inverse . ishift2D $ pad_mid a2 n
        a1a2 = zipWith (*) a1fft a2fft
        convolved = shift2D . myfft2D Forward $ a1a2
        mid = extract_mid n_ convolved
    in map (*n2) mid

-- More precise method
convolve2d' :: (Mode -> Acc (Matrix Visibility) -> Acc (Matrix Visibility))
            -> Acc (Matrix Visibility) -> Acc (Matrix Visibility) -> Acc (Matrix Visibility)
convolve2d' myfft a1 a2 =
    let
        (Z :. n :. _) = (unlift . shape) a1 :: Z :. Exp Int :. Exp Int
        m_ = 2*n - 1
        -- We need m to be a power of 2
        pw = A.ceiling (logBase 2 (A.fromIntegral m_) :: Exp F) :: Exp Int
        m = 2 ^ pw
        m2 = A.fromIntegral $ m * m

        f a = myfft Inverse . ishift2D $ pad_mid a m
        a1fft = f a1
        a2fft = f a2

        -- a1fft = myfft Inverse . ishift2D $ pad_mid a1 m
        -- a2fft = myfft Inverse . ishift2D $ pad_mid a2 m
        
        a1a2 = zipWith (*) a1fft a2fft
        convolved = shift2D . myfft Forward $ a1a2
        mid = extract_mid n convolved
    in map (*m2) mid

convolve2d :: Acc (Matrix Visibility) -> Acc (Matrix Visibility) -> Acc (Matrix Visibility)
convolve2d = convolve2d' myfft2D