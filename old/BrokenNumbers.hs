{-# language FlexibleContexts    #-}
{-# language RebindableSyntax    #-}
{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}

module BrokenNumbers where

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A
import Data.Array.Accelerate.Math.DFT.Centre              as A hiding (shift1D, shift2D, shift3D)
import Data.Array.Accelerate.Math.FFT                     as A

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..))
import Data.Char (isSpace)
import Text.Printf
import Data.List (intercalate)

import Debug.Trace
import System.IO.Unsafe (unsafePerformIO)
import Control.Exception (assert)

import Control.Lens as L (_1, _2, (^.) )

testFixbounds :: (Num a, Elt a) => Exp a -> Acc (Vector (Int, Int, a)) -> Acc (Matrix a)
testFixbounds cc0 v1 = 
    let 
        width = constant 5 :: Exp Int
        height = constant 5 :: Exp Int
        a = fill (index2 height width) cc0

        fullrepli origin = 
            let (x, y, source) = unzip3 $ v1 
                indexer id = 
                    let y' = y ! id
                        x' = x ! id
                    in lift (Z :. y' :. x')
            in permute (+) origin indexer source
    in fullrepli $ fullrepli a

c0 :: Exp Int
c0 = constant 0

vcomplex ::  Acc (Vector (Int, Int, Complex Double))
vcomplex = use $ fromList (Z :. 10) [((2 * x)  `P.mod` 5 ,(3 * x + 1) `P.mod` 5, ( (fromIntegral $ x+5) :+ 1.0)) | x <- [0..] ]

c0complex :: Exp (Complex Double)
c0complex =constant $  0.0 :+ 0.0

vdouble :: Acc (Vector (Int, Int, Double))
vdouble = use $ fromList (Z :. 10) [((2 * x)  `P.mod` 5 ,(3 * x + 1) `P.mod` 5, ( (fromIntegral $ x+5))) | x <- [0..] ]

c0double :: Exp Double
c0double = constant 0.0

test1 = I.run $ testFixbounds c0complex vcomplex
test1' = CPU.run $ testFixbounds c0complex vcomplex
test2 = I.run $ testFixbounds c0double vdouble
test2' = CPU.run $ testFixbounds c0double vdouble

main = do
    P.putStrLn "test1"
    P.putStrLn (P.show test1)
    P.putStrLn "test1'"
    P.putStrLn (P.show test1')
    P.putStrLn "test2"
    P.putStrLn (P.show test2)
    P.putStrLn "test2'"
    P.putStrLn (P.show test2')

{-
OUTPUT:
*BrokenNumbers> test1'
Matrix (Z :. 5 :. 5)
  [                                   0.0 :+ 0.0,                              42.0 :+ 4.0,                                 0.0 :+ 0.0,                                   0.0 :+ 0.0,                                    0.0 :+ 0.0,
                                     30.0 :+ 4.0,                               0.0 :+ 0.0,                                 0.0 :+ 0.0,                                   0.0 :+ 0.0,                                    0.0 :+ 0.0,
                                      0.0 :+ 0.0,                               0.0 :+ 0.0,                  0.0 :+ 1.40079706169e-312,  1.40079706181e-312 :+ 6.91046189490085e-310,                                   19.0 :+ 2.0,
     1.40079705984e-312 :+ 6.91046189490085e-310, 1.40079711786e-312 :+ 1.40079706189e-312, 1.400797059954e-312 :+ 6.910462099163e-310,                                  23.0 :+ 2.0, 6.91046646446026e-310 :+ 6.9104670154759e-310,
    6.91046189490085e-310 :+ 1.400797064183e-312, 1.40079706257e-312 :+ 1.40079706245e-312,                                17.0 :+ 2.0, 6.91046189490085e-310 :+ 1.400797063986e-312,  1.400797062637e-312 :+ 6.91046189490085e-310]
*BrokenNumbers> test1
Matrix (Z :. 5 :. 5)
  [  0.0 :+ 0.0, 42.0 :+ 4.0,  0.0 :+ 0.0,  0.0 :+ 0.0,  0.0 :+ 0.0,
    30.0 :+ 4.0,  0.0 :+ 0.0,  0.0 :+ 0.0,  0.0 :+ 0.0,  0.0 :+ 0.0,
     0.0 :+ 0.0,  0.0 :+ 0.0,  0.0 :+ 0.0,  0.0 :+ 0.0, 38.0 :+ 4.0,
     0.0 :+ 0.0,  0.0 :+ 0.0,  0.0 :+ 0.0, 46.0 :+ 4.0,  0.0 :+ 0.0,
     0.0 :+ 0.0,  0.0 :+ 0.0, 34.0 :+ 4.0,  0.0 :+ 0.0,  0.0 :+ 0.0]

*BrokenNumbers> test2'
Matrix (Z :. 5 :. 5)
  [  0.0, 42.0,  0.0,  0.0,  0.0,
    30.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0, 38.0,
     0.0,  0.0,  0.0, 46.0,  0.0,
     0.0,  0.0, 34.0,  0.0,  0.0]
*BrokenNumbers> test2
Matrix (Z :. 5 :. 5)
  [  0.0, 42.0,  0.0,  0.0,  0.0,
    30.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0, 38.0,
     0.0,  0.0,  0.0, 46.0,  0.0,
     0.0,  0.0, 34.0,  0.0,  0.0]
-}