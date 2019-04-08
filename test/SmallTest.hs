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
import Data.Array.Accelerate.Debug

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

processOne2 :: Acc (Scalar (Int, Int, Int)) -> Acc (Scalar(Int, Int, Int, Int, Visibility)) -> Acc (Array DIM5 Visibility) -> Acc (Array DIM3 Visibility) -> Acc (Matrix Visibility) -> Acc (Matrix Visibility) 
processOne2 _ _ _ _ t = t
-}
wkerns :: Acc (Array DIM5 Visibility)
wkerns = generate (constant (Z :. 1 :. 1 :. 1 :. 15 :. 15)) f
    where
        f :: Exp DIM5 -> Exp Visibility
        f (unlift -> (Z :. n :. yf :. xf :. y :. x) :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int :. Exp Int) = A.fromIntegral x * constant (0.01 :+ 0.005) + 0.1

akerns :: Acc (Array DIM3 Visibility)
akerns = generate (constant (Z :. 3 :. 15 :. 15)) f
    where
        f :: Exp DIM3 -> Exp Visibility
        f (unlift -> (Z :. n :. y :. x) :: Z :. Exp Int :. Exp Int :. Exp Int) = A.fromIntegral (x+y) * constant (0.01 :+ (-0.005)) + 0.008 + A.fromIntegral n * constant (0 :+ 0.1)

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

convTest2 :: Acc (Vector (Int, Int, Visibility))-- Acc (Matrix Visibility)--Acc (Vector (Int, Int, Int, Int))
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

        result1 :: Acc (Matrix Visibility)
        result1 = afor n processer2 (dest)
        

        -- The sequences stuff
        (id1, id2, id3) = A.unzip3 index
        visIndex = A.zip6 id1 id2 id3 cxf cyf vis
        
        visSeq :: Seq [Scalar (Int, Int, Int, Int, Int, Visibility)]
        visSeq = toSeq (constant (Z :. (0::Int))) visIndex

        coordsSeq :: Seq [Scalar (Int, Int)]
        coordsSeq = toSeq (constant (Z :. (0::Int))) $ A.zip cx cyf

        
        visKernSeq :: Seq [Matrix Visibility]
        visKernSeq = mapSeq (processOne2 wkerns akerns) visSeq
        {-
        visKernSeq :: Seq [Matrix Visibility]
        visKernSeq = mapSeq myProcessor visSeq
        -}
        myProcessor :: Acc (Scalar (Int, Int, Int, Int, Int, Visibility)) -> Acc (Matrix Visibility)
        myProcessor input = let
            gw = 15 
            gh = 15
            (_,_,_,_,_,vis) = unlift . the $ input :: (Exp Int, Exp Int, Exp Int, Exp Int, Exp Int, Exp Visibility)
            in fill (constant (Z :. gh :. gw)) vis
            

        addCoords :: Acc (Matrix Visibility) -> Acc (Scalar (Int, Int)) -> Acc (Matrix (Int, Int, Visibility))
        addCoords vis xy = let
            Z :. gh :. gw = unlift (shape vis) :: Z :. Exp Int :. Exp Int
            (x,y) = unlift . the $ xy :: (Exp Int, Exp Int)
            halfgw = gw `div` 2
            halfgh = gh `div` 2

            indexmapper :: Exp DIM2 -> Exp Visibility -> Exp (Int, Int, Visibility)
            indexmapper (unlift . unindex2 -> (i, j)::(Exp Int,Exp Int)) vis =
                    lift ( x + j - halfgw, y + i - halfgh, vis)

            in imap indexmapper vis

        visAndCoordSeq = zipWithSeq addCoords visKernSeq coordsSeq

        res :: Acc (Vector (Int, Int, Visibility))
        res = collect . elements $ visAndCoordSeq
        (xs, ys, resvis) = A.unzip3 res

        indexer id =
                let y' = ys ! id
                    x' = xs ! id
                in index2 y' x'

        result2 = permute (+) dest indexer resvis
    in res -- processer2 1 dest --result2

--res :: Acc (Scalar F)
--res = A.maximum . A.map real $ convTest2

res2 :: Acc (Scalar F)
res2 = A.maximum . A.map (\(unlift -> (x,y,v) :: (Exp Int, Exp Int, Exp Visibility)) -> real v) $ convTest2

ffttest :: Acc (Scalar Int) -> Acc (Scalar F)
ffttest i' = let i = the i'
    in A.maximum . A.map real . myfft2D Forward . (`pad_mid`32) $ slice akerns (lift (Z :. i :. All :. All))

ffttest2 = collect . elements $ produce 3 ffttest