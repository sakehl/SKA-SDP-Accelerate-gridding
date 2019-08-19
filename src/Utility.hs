{-# language FlexibleContexts    #-}
{-# language RebindableSyntax    #-}
{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}
module Utility where

import Types

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..), maybe)


pad_mid :: Acc (Matrix Visibility)    -- The input far field. Should be smaller than NxN
        -> Exp Int                    -- The desired far field size
        -> Acc (Matrix Visibility)    -- The far field, NxN
pad_mid ff n =
    let
        Z :. n0 :. n0w = (unlift . shape) ff :: Z :. Exp Int :. Exp Int
        result = if n == n0 then ff else padded
        pad_width = lift ((n `div` 2)-(n0 `div` 2), ((n+1) `div` 2)-((n0+1) `div` 2))
        padded = padder ff pad_width pad_width 0
    in result

--Extract a section from middle of a map Suitable for zero frequencies at N/2. This is the reverse operation to pad.
extract_mid :: Exp Int                  -- size of section
            -> Acc (Matrix Visibility)  -- grid from which to extract
            -> Acc (Matrix Visibility)
extract_mid n a = 
    let
        Z :. x :. y = (unlift. shape) a :: Z :. Exp Int :. Exp Int
        cx = x `div` 2
        cy = y `div` 2
        s  = n `div` 2
        --o  = if odd n then 1 else 0
        res = slit (cx -s) n . slit2 (cy -s) n $ a
        --r1 = dropOn _1 (cx -s) . takeOn _1 (cx + s + o) $ a
        --res = dropOn _2 (cy -s) . takeOn _2 (cy + s + o) $ r1
    in res  

extract_oversampled :: Acc (Matrix Visibility)        -- grid from which to extract
                    -> Exp Int                              -- oversampling factor
                    -> Exp Int                              -- size of section
                    -> Acc (Array DIM4 Visibility)    -- Return the oversampled kernels with shape qpx x qpx x n x n
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

------------------------
-- Helper functions
div3 :: Exp BaseLines -> Exp F -> Exp BaseLines
div3 (unlift -> (a,b,c)) x = lift (a / x, b / x, c / x)

myfor :: Int -> (Int -> a -> a) -> a -> a
myfor n f x | n P.== 0    = x
            | P.otherwise =  myfor (n-1) f (f (n-1) x)

liftTupf :: (Elt a, Elt b) => (Exp a -> Exp a) -> (Exp b -> Exp b) -> Exp (a,b) -> Exp (a,b)
liftTupf f g (unlift -> (a, b)) = lift (f a, g b)

afor :: forall a. Arrays a => Exp Int -> (Exp Int -> Acc a -> Acc a) -> Acc a -> Acc a
afor n f m = let
    newf :: Acc (Scalar Int, a) -> Acc (Scalar Int, a)
    newf (unlift -> (i, m) :: (Acc (Scalar Int), Acc a)) = 
        let i' = map (+1) i 
        in lift (i', f (the i) m)

    condition :: Acc (Scalar Int, a) -> Acc (Scalar Bool)
    condition (unlift -> (i, m) :: (Acc (Scalar Int), Acc a)) = map (<n) i

    initial :: Acc (Scalar Int, a)
    initial = lift (unit 0, m)
    in asnd $ awhile condition newf initial


padder :: Elt a => Acc (Matrix a) -> Exp (Int, Int) -> Exp (Int, Int) -> Exp a ->  Acc (Matrix a)
padder array pad_width_x pad_width_y constant_val =
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

c0 :: Exp Int
c0 = constant 0

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

findClosest :: (Elt a, Num a, Ord a) => Acc (Vector a) -> Exp a -> Exp Int
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



reverse1 :: Elt e => Acc (Matrix e) -> Acc (Matrix e)
reverse1 xs =
  let Z :. leny :. lenx = unlift . shape $ xs ::  Z :. Exp Int :. Exp Int
      pf id = let Z :. y :. x = unlift id :: Z :. Exp Int :. Exp Int
              in lift (Z :. y :. lenx - x -1)
  in  backpermute (shape xs) pf xs

reverse2 :: Elt e => Acc (Matrix e) -> Acc (Matrix e)
reverse2 xs =
  let Z :. leny :. lenx = unlift . shape $ xs ::  Z :. Exp Int :. Exp Int
      pf id = let Z :. y :. x = unlift id :: Z :. Exp Int :. Exp Int
              in lift (Z :. leny - y + 1 :. x)
  in  backpermute (shape xs) pf xs

slit2 :: Elt e => Exp Int -> Exp Int -> Acc (Matrix e) -> Acc (Matrix e)
slit2 m n acc =
  let m'        = the (unit m)
      n'        = the (unit n)
      Z :. sy :. sx  = unlift (shape acc)            :: Z :. Exp Int :. Exp Int
      index ix  = let Z :. j :. i = unlift ix        :: Z :. Exp Int :. Exp Int
                  in  lift (Z :. j + m' :. i)
  in
  backpermute (lift (Z :. n' `min` ((sy - m') `max` 0) :. sx ::Z :. Exp Int :. Exp Int)) index acc