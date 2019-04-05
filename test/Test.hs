{-# language FlexibleContexts    #-}
{-# language RebindableSyntax    #-}
{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A
import Data.Array.Accelerate.Math.DFT.Centre              as A
import Data.Array.Accelerate.Math.FFT                     as A

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I
import Data.Array.Accelerate.Debug                        as A

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..))

import Gridding

test :: Acc (Matrix Int) -> Acc (Array DIM3 Int, Matrix Int)
test pss = 
    let Z :. rows :. cols = unlift (shape pss) :: Z :. Exp Int :. Exp Int
        rss = fold (+) 0 . scanl (+) 0 $ pss
        rss' = replicate (lift (Z :. rows :. All)) rss
        rsss' = replicate (lift (Z :. rows :. rows :. All)) rss
        psss' = replicate (lift (Z :. rows :. All :. All)) pss
        asss = zipWith (+) psss' rsss'

        dss = fold (+) 0 asss
        wssN wss = zipWith (\w d -> 2*(d+w)) wss dss
        bss = wssN pss
        test = myfor 10 (\_ -> wssN) pss
    in lift (asss, bss)


tester :: Acc (Vector Double)
tester = use $ fromList (Z :. 2) [98..]

findClosest :: Acc (Vector Double) -> Exp Double -> Exp Int
findClosest ws w =
    let cmp (unlift -> (min, max) :: (Exp Int, Exp Int)) = (max - min) `div` 2 >= 1
        f (unlift -> (min, max) :: (Exp Int, Exp Int)) = 
            let id = (max + min) `div` 2
            in if w > ws !! id 
                then lift (id, max)
                else lift (min, id)
        Z :. max = unlift . shape $ ws :: Z :. Exp Int  
        minmax = lift (0 :: Exp Int, max)
        
        (r1, r2) = unlift . while cmp f $ minmax :: (Exp Int, Exp Int)
    in if abs (w - ws !! r1) < abs (w - ws !! r2) then r1 else r2