{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}
{-# language FlexibleContexts    #-}

module SmallTest where

import Types
import Gridding

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A
import Data.Array.Accelerate.Math.DFT.Centre              as A
import Data.Array.Accelerate.Math.FFT                     as A

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I

import Prelude as P
import Control.Lens as L (_1, _2, _3, _4, _5)


test :: Matrix Visibility
test = let processer = CPU.runN SmallTest.processOne2
           index0 i = CPU.run . unit $ index A.!! i
           dat0 i = CPU.run . unit $ dat A.!! i
           processi i = processer (index0 (constant i)) (dat0 (constant i)) wkerns' akerns'
           wkerns' = fromList (Z :. 0 :. 0 :. 0 :. 0 :. 0) []
           akerns' = fromList (Z :. 0 :. 0 :. 0) []
           a' = fromList (Z :. 1 :. 1) [1 :+ 1]

           index = use $ fromList (Z :. 4) [(x,x,x) | x <- [0..]] :: Acc (Vector (Int, Int, Int))
           dat = use $ fromList (Z :. 4) [(x,x,x,x, fromIntegral x) | x <- [0..]] :: Acc (Vector (Int, Int, Int, Int, Visibility))

       in (processi 3 P.$! processi 2 P.$! processi 1 P.$! processi 0 a')

processOne2 :: Acc (Scalar (Int, Int, Int)) -> Acc (Scalar(Int, Int, Int, Int, Visibility)) -> Acc (Array DIM5 Visibility) -> Acc (Array DIM3 Visibility) -> Acc (Matrix Visibility) -> Acc (Matrix Visibility) 
processOne2 _ _ _ _ t = t

{-
convgrid3 :: Runners (
             Acc (Array DIM5 Visibility)        -- The oversampled convolution w-kernel
          -> Acc (Array DIM3 Visibility)        -- The a-kernels
          -> Acc (Matrix Visibility)            -- Destination grid N x N of complex numbers
          -> Acc (Vector BaseLines)             -- The uvw baselines, but between -.5 and .5
          -> Acc (Vector (Int, Int, Int))       -- *DIF* The wbin index of the convolution kernel and the index of the a-kernels
          -> Acc (Vector Visibility)            -- The visiblities
          -> Acc (Matrix Visibility)
-}
wkerns :: Acc (Array DIM5 Visibility)
wkerns = generate (constant (Z :. 1 :. 1 :. 1 :. 15 :. 15)) f
    where
        f :: Exp DIM5 -> Exp Visibility
        f (unlift -> (Z :. n :. yf :. xf :. y :. x) :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int :. Exp Int) = A.fromIntegral x * constant (0.01 :+ 0.005) + 0.1

akerns :: Acc (Array DIM3 Visibility)
akerns = fill (constant (Z :. 3 :. 15 :. 15)) 0.1
    where
        f :: Exp DIM3 -> Exp Visibility
        f (unlift -> (Z :. n :. y :. x) :: Z :. Exp Int :. Exp Int :. Exp Int) = A.fromIntegral x * constant (0.01 :+ (-0.005)) + 0.008 + A.fromIntegral n * constant (0 :+ 0.1)

dest :: Acc (Matrix Visibility)
dest = fill (constant (Z :. 10 :. 10)) 0

uvw :: Acc (Vector BaseLines)
uvw = use $ fromList (Z :. 2) [(0.1, 0.2, 0.3), (-0.1, 0.4, 0.1)]

index :: Acc (Vector (Int, Int, Int))
index = use $ fromList (Z :. 2) [(0, 0, 1), (0, 0, 2)]

vis :: Acc (Vector Visibility)
vis = use $ fromList (Z :. 2) [0.3 :+ 0.5, 0.4 :+ 0.2]

convTest :: Acc (Matrix Visibility)
convTest = convgrid3 wkerns akerns dest uvw index vis

convTest2 ::  Acc (Matrix Visibility)--Acc (Vector (Int, Int, Int, Int))
convTest2 =
    let Z :. _ :. qpx :. _ :. gh :. gw = unlift (shape wkerns) :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int :. Exp Int
        Z :. height :. width = unlift (shape dest) :: Z :. Exp Int :. Exp Int
        Z :. n = unlift (shape vis) :: Z :. Exp Int
        halfgh = gh `div` 2
        halfgw = gw `div` 2

        --Gather the coords, make them into ints and fractions and zip with visibility
        coords = frac_coords (lift (height, width)) qpx uvw
        (cx, cxf, cy, cyf) = unzip4 coords
        dat = zip5 cx cxf cy cyf vis

        index03 = index
        dat03   = dat
        wkerns03 = wkerns
        akerns03 = akerns
        processer2 :: Exp Int -> Acc (Matrix Visibility) -> Acc (Matrix Visibility)
        processer2 n aa = 
            let id = unit $ index03 A.!! n
                dat = unit $ dat03 A.!! n
            in processOne id dat wkerns03 akerns03 aa

        result2 :: Acc (Matrix Visibility)
        result2 = afor n processer2 (dest)
    in result2

convTest3 ::  Acc (Matrix Visibility)--Acc (Vector (Int, Int, Int, Int))
convTest3 =
    let Z :. _ :. qpx :. _ :. gh :. gw = unlift (shape wkerns) :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int :. Exp Int
        Z :. height :. width = unlift (shape dest) :: Z :. Exp Int :. Exp Int
        Z :. n = unlift (shape vis) :: Z :. Exp Int
        halfgh = gh `div` 2
        halfgw = gw `div` 2

        --Gather the coords, make them into ints and fractions and zip with visibility
        coords = frac_coords (lift (height, width)) qpx uvw
        (cx, cxf, cy, cyf) = unzip4 coords
        (wbin, a1, a2) = A.unzip3 index
        visIndex = A.zip6 wbin a1 a2 cxf cyf vis

        addCoords :: Acc (Matrix Visibility) -> Acc (Scalar (Int, Int)) -> Acc (Matrix (Int, Int, Visibility))
        addCoords vis xy = let
            (x,y) = unlift . the $ xy :: (Exp Int, Exp Int)
            indexmapper (unlift . unindex2 -> (i, j)::(Exp Int,Exp Int)) vis =
                    lift ( x + j - halfgw, y + i - halfgh, vis)
            in imap indexmapper vis
        
        empty :: Acc (Array DIM3 Visibility)
        empty = fill (lift $ Z :. (0::Exp Int) :. gh :. gw) 0

        processer :: Exp Int -> Acc (Array DIM3 Visibility) -> Acc (Array DIM3 Visibility)
        processer n aa = 
            let id = unit $ index A.!! n
                visId = unit $ visIndex A.!! n
                res = processOne3 wkerns akerns visId
            in concatMat aa res 

        kernVis :: Acc (Array DIM3 Visibility)
        kernVis = afor n processer empty
        
        createAllCoords :: Acc (Vector (Int, Int)) -> Acc (Array DIM3 (Int, Int))
        createAllCoords xy = let
            replxy = A.replicate (lift $ Z :. All :. gh :. gw) xy
            mapper (unlift -> Z :. n :. j :. i :: Z :. Exp Int :. Exp Int :. Exp Int) (unlift -> (x,y) :: (Exp Int, Exp Int)) =
                lift ( x + i - halfgw, y + j - halfgh)
            in imap mapper replxy
        
        (x',y') = A.unzip . createAllCoords $ A.zip cx cy
        withCoords = flatten $ A.zip3 x' y' kernVis

        --Fix the bounds
        fixer = fixoutofbounds width height 0
        fixedBounds = A.map fixer withCoords
        (xs,ys, val) = A.unzip3 fixedBounds

        indexer id =
            let y' = ys ! id
                x' = xs ! id
            in index2 y' x'
    in permute (+) dest indexer val


processOne3 :: Acc (Array DIM5 Visibility) -> Acc (Array DIM3 Visibility) -> Acc (Scalar (Int, Int, Int, Int, Int, Visibility)) -> Acc (Matrix Visibility)
processOne3 wkerns akerns 
    (unlift . the -> (wbin, a1index, a2index, xf, yf, vis) :: (Exp Int, Exp Int, Exp Int, Exp Int, Exp Int, Exp Visibility)) =
        let Z :. _ :. _ :. _ :. gh :. gw = unlift (shape wkerns) :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int :. Exp Int
            halfgh = gh `div` 2
            halfgw = gw `div` 2
            
            --Get the right a and w kernels
            a1 = slice akerns (lift (Z :. a1index :. All :. All))
            a2 = slice akerns (lift (Z :. a2index :. All :. All))
            w  = slice wkerns (lift (Z :. wbin :. All :. All :. All :. All))
            --Convolve them
            -- NOTE, the conjugate normally happens in imaging function, but it is convenient to do it here.
            awkern = A.map conjugate $ aw_kernel_fn2 yf xf w a1 a2

            -- Let the visibility have the same dimensions as the aw-kernel
            allvis = A.replicate ( lift ( Z :. gh :. gw)) (unit vis)
        in A.zipWith (*) allvis awkern