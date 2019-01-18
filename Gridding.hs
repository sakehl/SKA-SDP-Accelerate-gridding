{-# language FlexibleContexts    #-}
{-# language RebindableSyntax    #-}
{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}
-- {-# language ScopedTypeVariables #-}

module Gridding where

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), IO)

-------------------------------
-- Gridding Accelerate code

simple_imaging :: Exp Float                             -- Field of view size
               -> Exp Int                               -- Grid size
               -> Vector (Float, Float, Float)          -- all the uvw baselines (coordinates) (lenght : n * (n-1))
               -> Vector (Int, Int, Float, Float)       -- NOT USED HERE (Antenna 1, Antenna 2, The time (in MJD UTC), Frequency (Hz)) (lenght : n * (n-1))
               -> Vector (Complex Float)                -- visibility  (length n * (n-1))
               -> Acc (Matrix (Complex Float))
simple_imaging theta lam uvw src vis =
    let n = round(theta * fromIntegral lam)
        p = map (`div3` fromIntegral lam) (use uvw)
        czero = constant $ 0 :+ 0
        guv = fill (index2 n n) czero
    in grid guv p (use vis)
    where
        div3 :: Exp (Float, Float, Float) -> Exp Float -> Exp (Float, Float, Float)
        div3 (unlift -> (a,b,c)) x = lift (a / x, b / x, c / x)

grid :: Acc (Matrix (Complex Float))                    --Destination grid N x N of complex numbers
     -> Acc (Vector (Float, Float, Float))              -- The uvw baselines, but between -.5 and .5
     -> Acc (Vector (Complex Float))                    -- The visiblities
     -> Acc (Matrix (Complex Float))
grid a p v = permute (+) a (\i -> xy ! i) v
    where
        Z :. n :. _ = unlift (shape a) :: Z :. Exp Int :. Exp Int
        halfn = n `div` 2
        nf = fromIntegral n :: Exp Float
        xy = map gridxy p :: Acc (Vector DIM2)

        gridxy :: Exp (Float, Float, Float) -> Exp DIM2
        gridxy (unlift -> (x, y, _ :: Exp Float)) = index2 (toGridCell x) (toGridCell y)

        toGridCell :: Exp Float -> Exp Int
        toGridCell f = halfn + floor (0.5 + nf * f)


-------------------------------
-- Gridding Python code
{-
def simple_imaging(--Field of view size
                    theta :: Float
                     --Grid size
                  , lam  :: Int
                  -- all the uvw baselines (coordinates) (lenght : n * (n-1))                          
                  , uvw  :: [(Float, Float, Float)]
                  -- (Antenna 1, Antenna 2, The time (in MJD UTC), Frequency (Hz)) (lenght : n * (n-1))                     
                  , src  :: [(Int, Int, Float, Float)] 
                  -- The visibility 
                  , vis :: [Complex Float]
                  ):
    """Trivial function for imaging
    Does no convolution but simply puts the visibilities into a grid
    cell i.e. boxcar gridding
    """
    N = int(round(theta * lam))
    assert N > 1
    //The uv plane to grid to N x N matrix of complex numbers
    guv = numpy.zeros([N, N], dtype=complex)
    grid(guv, uvw / lam, vis)
return guv


def grid(
        --Destination grid N x N of complex numbers
         a :: Matrix (Complex Float)
         -- The uvw baselines, but between -.5 and .5
         , p :: [(Float, Float, Float)]
        -- The visiblities
         , v :: [Complex Float]
         ):
    """Grid visibilities (v) at positions (p) into (a) without convolution
    :param a:   The uv plane to grid to (updated in-place!)
    :param p:   The coordinates to grid to (in fraction [-.5,.5[ of grid)
    :param v:   Visibilities to grid
    """
    assert numpy.max(p) < 0.5

    N = a.shape[0]
    xy = N//2 + numpy.floor(0.5 + N * p[:,0:2]).astype(int)
    for (x, y), v in zip(xy, v):
        a[y, x] += v
-}