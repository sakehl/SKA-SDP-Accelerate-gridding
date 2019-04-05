{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}

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

{-
test :: Matrix Visibility
test = let processer = undefined --CPU.runN SmallTest.processOne2
           index0 i = CPU.run . unit $ index A.!! i
           dat0 i = CPU.run . unit $ dat A.!! i
           processi i = processer (index0 (constant i)) (dat0 (constant i)) wkerns' akerns'
           wkerns' = fromList (Z :. 0 :. 0 :. 0 :. 0 :. 0) []
           akerns' = fromList (Z :. 0 :. 0 :. 0) []
           a' = fromList (Z :. 1 :. 1) [1 :+ 1]

           index = use $ fromList (Z :. 4) [(x,x,x) | x <- [0..]] :: Acc (Vector (Int, Int, Int))
           dat = use $ fromList (Z :. 4) [(x,x,x,x, fromIntegral x) | x <- [0..]] :: Acc (Vector (Int, Int, Int, Int, Visibility))

       in (processi 3 P.$! processi 2 P.$! processi 1 P.$! processi 0 a')
-}
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

    in processer2 0 dest --result2
{-
convgrid3 :: Runners (
             Acc (Array DIM5 Visibility)        -- The oversampled convolution w-kernel
          -> Acc (Array DIM3 Visibility)        -- The a-kernels
          -> Acc (Matrix Visibility)            -- Destination grid N x N of complex numbers
          -> Acc (Vector BaseLines)             -- The uvw baselines, but between -.5 and .5
          -> Acc (Vector (Int, Int, Int))       -- *DIF* The wbin index of the convolution kernel and the index of the a-kernels
          -> Acc (Vector Visibility)            -- The visiblities
          -> Acc (Matrix Visibility)
        )
convgrid3 run runN wkerns akerns a p index v =
    let Z :. _ :. qpx :. _ :. gh :. gw = unlift (shape wkerns) :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int :. Exp Int
        Z :. height :. width = unlift (shape a) :: Z :. Exp Int :. Exp Int
        Z :. n = unlift (shape v) :: Z :. Exp Int
        halfgh = gh `div` 2
        halfgw = gw `div` 2

        --Gather the coords, make them into ints and fractions and zip with visibility
        coords = frac_coords (lift (height, width)) qpx p
        (cx, cxf, cy, cyf) = unzip4 coords
        dat = zip5 cx cxf cy cyf v

        -- We reuse this one a lot, so compile it
        {-
        processer = runN processOne
        index02 = runN (\i -> unit $ index !! the i)
        dat02   = runN (\i -> unit $ dat   !! the i)

        dat0 i = dat02 (fromList Z [i])
        index0 i = index02 (fromList Z [i])
        wkerns' = run wkerns
        akerns' = run akerns
        a' = run a
        n' = P.head . toList .run . unit $  n
        {-
        processi i | i `P.mod` 1000 P.== 0 = trace (printf "We are processing number %i now (total: %i)" t0) i n') 
                            $ processer (index0 i) (dat0 i) wkerns' akerns'
                   | otherwise = processer (index0 i) (dat0 i) wkerns' akerns'-}
        processi i = processer (index0 i) (dat0 i) wkerns' akerns'

        result = myfor n' processi a'
        -}

        index03 = index
        dat03   = dat
        wkerns03 = wkerns
        akerns03 = akerns
        processer2 :: Exp Int -> Acc (Matrix Visibility) -> Acc (Matrix Visibility)
        processer2 n aa = 
            let id = unit $ index03 !! n
                dat = unit $ dat03 !! n
            in processOne id dat wkerns03 akerns03 aa

        result2 :: Acc (Matrix Visibility)
        result2 = afor n processer2 (a)

        
        
    in result2

-}