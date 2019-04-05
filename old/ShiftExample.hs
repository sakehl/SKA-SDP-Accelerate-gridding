{-# language FlexibleContexts    #-}
{-# language RebindableSyntax    #-}
{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}


module ShiftExample where

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I

import Data.Array.Accelerate.Math.DFT.Centre              as A

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO)
import Data.Char
import Text.Printf
import Data.List (intersperse, intercalate)

import Debug.Trace
import System.IO.Unsafe (unsafePerformIO)
import Control.Exception

ffttest0 :: Int -> (Vector (Double, Int) , Vector Double) -- , [(Double, Int)], Vector Int)
ffttest0 n = let
    halfn = fromIntegral $ (n - 1) `div` 2
    minstart = case P.even n of
        True -> -halfn - 1
        False -> -halfn
    lst = [0..halfn] P.++ [minstart..(-1)]
    freqs = fromList (Z :. n) lst

    res2 = shift1D' (use freqs)
    
    res4 = shiftHELP (use freqs)
    res5 = zip (use freqs) res4
    Z:.w    = unlift (A.shape (use freqs)) :: Z :. Exp Int
    in (CPU.run res5, CPU.run res2) --), res3, CPU.run res4)

shift1D' :: Elt e => Acc (Vector e) -> Acc (Vector e)
shift1D' arr
  = A.backpermute (A.shape arr) p arr
  where
    p ix
      = index1 $ theFunc arr ix

shiftHELP :: Elt e => Acc (Vector e) -> Acc (Vector Int)
shiftHELP arr
  = imap p arr
  where
    p ix _
      -- = indexHead $ (ilift1 P.id) ix
       = theFunc arr ix

theFunc3 :: Elt e => Acc (Vector e) -> Exp (Z :. Int) -> Exp Int
theFunc3 _ ix = indexHead $ (ilift1 P.id) ix

theFunc :: Elt e => Acc (Vector e) -> Exp (Z :. Int) -> Exp Int
theFunc arr ix =
    let Z:.x = unlift ix :: Z :. Exp Int
    in (x A.< mw ? (x + mw2, x - mw))
    where
        Z:.w    = unlift (shape arr) :: Z :. Exp Int
        mw      = w `div` 2
        mw2     = (w + 1) `div` 2

shiftTest :: [Double] -> [(Double, Int)]
shiftTest xs = let
    temp = addindex 0 xs
    addindex n [] = []
    addindex n (x:xs) = (x, n) : addindex (n+1) xs

    w = P.length xs
    mw      = (w) `P.div` 2
    mw2     = (w) `P.div` 2
    indexer (v,x) = case x P.< mw2 of
        True  -> (v, x + mw)
        False -> (v, x -mw2)
    
    boolToInt b = case b of
        False -> 0
        True -> 1
    
    n = P.length xs
    shift = (n `P.quot` 2) + boolToInt (P.odd n)
    roll i  = (i + shift) `P.rem` n
    
    indexer2 (v,x) = (v, roll x)
    in P.map indexer temp



shift1D :: Elt e => Acc (Vector e) -> Acc (Vector e)
shift1D arr = backpermute sh p arr
  where
    sh      = shape arr
    n       = indexHead sh
    --
    shift   = (n `quot` 2) + boolToInt (odd n)
    roll i  = (i+shift) `rem` n
    p       = ilift1 roll

theFunc2 :: Elt e => Acc (Vector e) -> Exp (Z :. Int) -> Exp Int
theFunc2 arr ix =
    let Z:.x = unlift ix :: Z :. Exp Int
    in roll x
    where
        sh      = shape arr
        n       = indexHead sh
        --
        shift   = (n `quot` 2) + boolToInt (odd n)
        roll i  = (i+shift) `rem` n
        p       = ilift1 roll