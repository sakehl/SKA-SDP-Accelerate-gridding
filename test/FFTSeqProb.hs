{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}
{-# language FlexibleContexts    #-}
{-# language RebindableSyntax    #-}
module FFTSeqProb where

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A
import Data.Array.Accelerate.Math.DFT.Centre              as A
import Data.Array.Accelerate.Math.FFT                     as A

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I
import Data.Array.Accelerate.Debug

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..), maybe)

type Visibility = Complex Double
type F = Double
type (Matrix a) = Array DIM2 a

akerns :: Acc (Array DIM3 Visibility)
akerns = generate (constant (Z :. 3 :. 16 :. 16)) f
    where
        f :: Exp DIM3 -> Exp Visibility
        f (unlift -> (Z :. n :. y :. x) :: Z :. Exp Int :. Exp Int :. Exp Int) = A.fromIntegral (x+y) * constant (0.01 :+ (-0.005)) + 0.008 + A.fromIntegral (n+1) * constant (0.2 :+ 0.1)

testMat :: Acc (Array DIM3 F)
testMat = generate (constant (Z :. 3 :. 16 :. 16)) f
    where
        f :: Exp DIM3 -> Exp F
        f (unlift -> (Z :. n :. y :. x) :: Z :. Exp Int :. Exp Int :. Exp Int) = A.fromIntegral (x+y) * constant 0.01 + 0.008 + A.fromIntegral (n+1) * constant 0.2

myfft2D :: FFTElt e => A.Mode -> Acc (Array DIM2 (Complex e)) -> Acc (Array DIM2 (Complex e))
myfft2D mode arr = let
    --sh = P.head . toList . CPU.run . unit$ shape arr
    sh = Z :. 32 :. 32
    in fft2D' mode sh arr

ffttest :: Acc (Scalar Int) -> Acc (Scalar F)
ffttest i' = let i = the i'
    in A.maximum . A.map real . myfft2D Forward . (`pad_mid`32) $ slice akerns (lift (Z :. i :. All :. All))

ffttest2 = collect . elements $ produce 3 ffttest


permuteTest :: Acc (Scalar Int) -> Acc (Scalar F)
permuteTest i' = let i = the i'
    in A.maximum . (`pad_mid`32) $ slice testMat (lift (Z :. i :. All :. All))

permuteTest2 = collect . elements $ produce 3 permuteTest

testF :: Acc (Scalar Int) -> Acc (Array DIM2 Int)
testF i' = let
    i = the i'
    orig = generate (index2 4 4) shapeSize
    def = fill (index2 5 5) i
 in permute const def P.id orig

permuteTest3 = collect . elements $ produce 3 testF


padder :: Elt a => Acc (Matrix a) -> Exp (Int, Int) -> Exp (Int, Int) -> Exp a ->  Acc (Matrix a)
padder array pad_width_x pad_width_y constant_val =
    let
        (x0, x1) = unlift pad_width_x :: (Exp Int, Exp Int)
        (y0, y1) = unlift pad_width_y :: (Exp Int, Exp Int)
        Z :. m :. n = (unlift . shape) array :: ( Z :. Exp Int :. Exp Int)
        def = fill (index2 (m + y0 + y1) (n + x0 + x1)) constant_val

        indexer (unlift -> Z :. y :. x :: Z :. Exp Int :. Exp Int)
            = index2 (y + y0) (x + x0)
        result = permute const def indexer array
    in result

padder2 :: Elt a => Acc (Matrix a) -> Exp (Int, Int) -> Exp (Int, Int) -> Exp a ->  Acc (Matrix a)
padder2 array pad_width_x pad_width_y constant_val =
    let
        (x0, x1) = unlift pad_width_x :: (Exp Int, Exp Int)
        (y0, y1) = unlift pad_width_y :: (Exp Int, Exp Int)
        Z :. m :. n = (unlift . shape) array :: ( Z :. Exp Int :. Exp Int)
        newshape = index2 (m + y0 + y1) (n + x0 + x1)

        indexer (unlift -> Z :. y :. x :: Z :. Exp Int :. Exp Int)
            = let oldx = x - x0
                  oldy = y - y0
                  inrange = oldx >= 0 && oldx < n && oldy >= 0 && oldy < m
                  oldval  = array ! index2 oldx oldy
            in if inrange then oldval else constant_val
    in generate newshape indexer

pad_mid :: (Elt a, Num a)=> Acc (Matrix a)    -- The input far field. Should be smaller than NxN
        -> Exp Int                    -- The desired far field size
        -> Acc (Matrix a)    -- The far field, NxN
pad_mid ff n =
    let
        Z :. n0 :. n0w = (unlift . shape) ff :: Z :. Exp Int :. Exp Int
        result = if n == n0 then ff else padded
        pad_width = lift ((n `div` 2)-(n0 `div` 2), ((n+1) `div` 2)-((n0+1) `div` 2))
        padded = padder2 ff pad_width pad_width 0
    in result