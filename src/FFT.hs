{-# LANGUAGE CPP                   #-}
{-# LANGUAGE ConstraintKinds       #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE GADTs                 #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE RankNTypes            #-}
{-# LANGUAGE RebindableSyntax      #-}
{-# LANGUAGE ScopedTypeVariables   #-}
{-# LANGUAGE TypeOperators         #-}
{-# LANGUAGE ViewPatterns          #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.Adhoc
-- Copyright   : [2017] Henning Thielemann
--               [2017] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- Implementation of ad-hoc FFT stolen from the accelerate-fourier by Henning
-- Thielemann (BSD3 licensed), and updated to work with current Accelerate. That
-- package contains other more sophisticated algorithms as well.
--

module FFT (myfft2D, ifft, fft, ditSplitRadixLoop, Numeric(..), NumericR(..))
  where

import           Data.Array.Accelerate                      as A hiding
                                                                  (transpose)
import           Data.Array.Accelerate.Array.Sugar          (shapeToList,
                                                             showShape)
import           Data.Array.Accelerate.Control.Lens.Shape
import           Data.Array.Accelerate.Data.Bits            as A
import           Data.Array.Accelerate.Data.Complex
import           Data.Array.Accelerate.Math.FFT             hiding
                                                             (ditSplitRadixLoop,
                                                             fft)
import qualified Data.Array.Accelerate.Math.FFT.LLVM.Native as Native
import qualified Data.Array.Accelerate.Math.FFT.LLVM.PTX    as PTX
import           Data.Bits                                  as B
import qualified Prelude                                    as P
import           Text.Printf


data NumericR a where
  NumericRfloat32 :: NumericR Float
  NumericRfloat64 :: NumericR Double

class (RealFloat a, FromIntegral Int a, Elt (Complex a)) => Numeric a where
  numericR :: NumericR a

instance Numeric Float where
  numericR = NumericRfloat32

instance Numeric Double where
  numericR = NumericRfloat64

{-
-- Somehow GHCI doesn't pick up Data.Array.Accelerate.Control.Lens.Shape, so we just redifine the lens instances here
import Control.Lens (lens)

instance (Slice sh, Elt a, Elt a')
    => Field1 (Exp (sh :. a)) (Exp (sh :. a')) (Exp a) (Exp a') where
  _1 = lens indexHead (\sh b -> lift (indexTail sh :. b))

instance (Slice sh, Elt a, Elt b, Elt b', Slice (sh :. b), Slice (sh :. b'))
    => Field2 (Exp (sh :. b :. a)) (Exp (sh :. b' :. a)) (Exp b) (Exp b') where
  _2 = lens (\s   -> let _  :. b :. _ = unlift s :: Exp sh :. Exp b :. Exp a in b)
            (\s b -> let sh :. _ :. a = unlift s :: Exp sh :. Exp b :. Exp a in lift (sh :. b :. a))
-}
myfft2D :: (Numeric e, P.Num e, IsFloating e)
    => Mode
    -> Acc (Array DIM2 (Complex e))
    -> Acc (Array DIM2 (Complex e))
myfft2D = myfft2DV3

myfft2DV1 :: (Numeric e, P.Num e, IsFloating e)
    => Mode
    -> Acc (Array DIM2 (Complex e))
    -> Acc (Array DIM2 (Complex e))
myfft2DV1 mode arr =
  let
    scale = A.fromIntegral (A.size arr)
    go    = ditSplitRadixLoop mode
  in case mode of
      Inverse -> A.map (/scale) (go arr)
      _       -> go arr

myfft2DV2 :: (Numeric e, P.Num e, IsFloating e)
    => Mode
    -> Acc (Array DIM2 (Complex e))
    -> Acc (Array DIM2 (Complex e))
myfft2DV2 mode arr =
  let
    scale  = A.fromIntegral (A.size arr)
    go     = transpose . fftV2 sign (Z:.32) width . transpose . fftV2 sign (Z:.width)  height
    height = 32
    width  = 32
    sign   = signOfMode mode
  in case mode of
      Inverse -> A.map (/scale) (go arr)
      _       -> go arr

myfft2DV3 :: (Numeric e, P.Num e, IsFloating e)
    => Mode
    -> Acc (Array DIM2 (Complex e))
    -> Acc (Array DIM2 (Complex e))
myfft2DV3 = fft2DFor

ifft :: (Shape sh, Slice sh, Numeric e)
    => Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
ifft = fft Inverse

fft :: (Shape sh, Slice sh, Numeric e)
    => Mode
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
fft mode arr =
  let len             = indexHead (shape arr)
      (pow2, smooth5) = is2or5smooth len
  in
  if len <= 1 then arr                        else
  if pow2     then ditSplitRadixLoop mode arr else
  if smooth5  then dit235            mode arr
              else transformChirp235 mode arr


-- Implementations
-- ---------------

is2or5smooth :: Exp Int -> (Exp Bool, Exp Bool)
is2or5smooth len =
  let maxPowerOfTwo = len A..&. negate len
      lenOdd        = len `quot` maxPowerOfTwo
  in
  ( 1 == lenOdd
  , 1 == divideMaxPower 5 (divideMaxPower 3 lenOdd)
  )

divideMaxPower :: Exp Int -> Exp Int -> Exp Int
divideMaxPower fac =
  while (\n -> n `rem`  fac == 0)
        (\n -> n `quot` fac)


-- -- | Split-radix for power-of-two sizes
-- --
-- ditSplitRadix
--     :: (Shape sh, Numeric e)
--     => Mode
--     -> Acc (Array (sh:.Int) (Complex e))
--     -> Acc (Array (sh:.Int) (Complex e))
-- ditSplitRadix mode arr =
--   if indexHead (shape arr) <= 1
--     then arr
--     else ditSplitRadixLoop mode arr

ditSplitRadixLoop
    :: forall sh e. (Shape sh, Slice sh, Numeric e)
    => Mode
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
ditSplitRadixLoop mode arr =
  let
      twiddleSR (fromIntegral -> n4) k (fromIntegral -> j) =
        let w = pi * k * j / (2*n4)
        in  lift (cos w :+ signOfMode mode * sin w)

      --Twidles with same len4, generate same shaped arrays
      twiddle len4 k =
        generate (index1 len4) (twiddleSR len4 k . indexHead)
      -- So step stays regular
      step (unlift -> (us,zs)) =
        let
            --k is same, since dependent on shape
            k           = indexHead (shape zs)
            --Regulary shaped
            tw1         = twiddle k 1
            --Regulary shaped
            tw3         = twiddle k 3
            --
            im          = lift (makeExp 0 :+ signOfMode mode)
            twidZeven   = zipWithExtrude1 (*) tw1 (sieveV 2 0 zs) --sieveV gives regular shape
            twidZodd    = zipWithExtrude1 (*) tw3 (sieveV 2 1 zs) -- since input is reguar, zipwithExtrude1 is aswell
            zsum        = zipWith (+) twidZeven twidZodd          -- input regular, thus also output
            zdiff       = map (im *) (zipWith (-) twidZeven twidZodd) --input regular, thus output
            zcomplete   = zsum ++ zdiff -- input regular, thus output
            _ :. n :. _ = unlift (shape zcomplete) :: Exp sh :. Exp Int :. Exp Int -- same for each, since regular
        in
        lift ( zipWith (+) us zcomplete ++ zipWith (-) us zcomplete -- This is regular
             , dropV n us
             )

      --transform stays regular, so everything stays regular
      rebase s = lift (transform2 (-1) (afst s), asnd s)

      --Reorder stays regular, for reasons below
      reorder (unlift -> (xs,ys)) =
        let -- Evens and odd have stay regular, since the factor (2) is always the same.
            evens = sieve 2 0 xs
            odds  = sieve 2 1 xs
        in
        -- Since evens stay regular and twist stays regular (factor is always 2), this stays regular
        lift (evens ++^ ys, twist 2 odds)
      -- New shape is comming from old shape, fully determined by it. So stays regular.
      initial :: Acc (Array (sh :. Int :. Int) (Complex e), Array (sh :. Int :. Int) (Complex e))
      initial =
        let sh :. n = unlift (shape arr) :: Exp sh :. Exp Int
        in  lift ( reshape (lift (sh :. (1::Int) :. n)) arr
                 , fill    (lift (sh :. (0::Int) :. n `quot` 2)) 0
                 )
  in
  --headV also stays regular, so result is regular
  headV
    $ afst
    -- step stay regular, and since check function is fully determined on the shape, it will stay regular
    $ awhile (\s -> unit (indexHead (indexTail (shape (asnd s))) > 0)) step
    -- rebase stays regular
    $ rebase
    -- Since the check function is fully determined on the shape, it will stay regular
    $ awhile (\s -> unit (indexHead (shape (asnd s)) > 1)) reorder
    $ initial


-- | Decimation in time for sizes that are composites of the factors 2,3 and 5.
-- These sizes are known as 5-smooth numbers or the Hamming sequence.
--
-- <http://oeis.org/A051037>
--
dit235
    :: forall sh e. (Shape sh, Slice sh, Numeric e)
    => Mode
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
dit235 mode arr =
  let
      merge :: forall sh' a. (Shape sh', Slice sh', Elt a)
            => Acc (Array (sh':.Int:.Int) a)
            -> Acc (Array (sh':.Int) a)
      merge xs =
        let sh :. m :. n = unlift (shape xs) :: Exp sh' :. Exp Int :. Exp Int
        in  backpermute
              (lift (sh :. m*n))
              (\(unlift -> ix :. k :: Exp sh' :. Exp Int) ->
                  let (q,r) = k `quotRem` m
                  in  lift (ix :. r :. q))
              xs

      step fac xs =
        let sh :. count :. len = unlift (shape xs) :: Exp sh :. Exp Int :. Exp Int
            --
            twiddled :: Acc (Array (sh:.Int:.Int:.Int) (Complex e))
            twiddled = transpose
                     $ zipWithExtrude2 (*) (twiddleFactors fac len)
                     $ reshape (lift (sh :. count `quot` fac :. fac :. len)) xs
        in
        merge $ if fac == 5 then transform5 cache5 twiddled else
                if fac == 4 then transform4 cache4 twiddled else
                if fac == 3 then transform3 cache3 twiddled
                            else transform2 cache2 twiddled

      initial :: Acc (Array (sh:.Int:.Int) (Complex e), Vector Int)
      initial =
        let sh :. n = unlift (shape arr) :: Exp sh :. Exp Int
        in  lift ( reshape (lift (sh :. (1::Int) :. n)) arr
                 , fill (index1 0) 0
                 )

      twiddleFactors :: Exp Int -> Exp Int -> Acc (Matrix (Complex e))
      twiddleFactors m n =
        generate (index2 m n)
                 (\(unlift -> Z :. j :. i) -> twiddle (m*n) j i)

      cisrat :: Exp Int -> Exp Int -> Exp (Complex e)
      cisrat d n =
        let w = 2*pi * fromIntegral n / fromIntegral d
        in  lift (cos w :+ signOfMode mode * sin w)

      twiddle :: Exp Int -> Exp Int -> Exp Int -> Exp (Complex e)
      twiddle n k j = cisrat n ((k*j) `rem` n)

      cache2 :: Exp (Complex e)
      cache2 = -1

      cache3 :: Exp (Complex e, Complex e)
      cache3 =
        let sqrt3d2 = sqrt 3 / 2
            mhalf   = makeExp (-1) / makeExp 2 :: Exp e-- lift (-1/2) :: Exp e
            s       = signOfMode mode
            u       = s * sqrt3d2
        in
        lift ((mhalf :+ u, mhalf :+ (-u)))

      cache4 :: Exp (Complex e, Complex e, Complex e)
      cache4 =
        let s = signOfMode mode :: Exp e

            --t1 = constant ((-1) :: e) :: Exp e
            t1 = makeExp (-1) :: Exp e
            t2 = makeExp 0 :: Exp e

            r1 = lift (t2 :+ s) :: Exp (Complex e)
            r2 = lift (t1 :+ t2) :: Exp (Complex e)
            r3 = lift (t2 :+ (-s)) :: Exp (Complex e)
        in lift (r1, r2, r3)
        --in lift (constant 0 :+ s, constant ((-1) :+ (-0)), constant 0 :+ (-s))
        --in lift (0 :+ s, (-1) :+ (-0), 0 :+ (-s))

      cache5 :: Exp (Complex e, Complex e, Complex e, Complex e)
      cache5 =
        let z = cisrat 5
        in  lift (z 1, z 2, z 3, z 4)
  in
  headV
    $ afst
    $ awhile
        (\s -> unit (length (asnd s) > 0))
        (\s -> let (xs,fs) = unlift s
                   f       = fs !! 0
               in
               lift (step f xs, tail fs))
    $ awhile
        (\s -> unit (indexHead (shape (afst s)) > 1))
        (\s -> let (xs,fs)      = unlift s
                   len          = indexHead (shape xs)
                   divides k n  = n `rem` k == 0
                   f            = if divides 3 len then 3 else
                                  if divides 4 len then 4 else
                                  if divides 5 len then 5
                                                   else 2
               in
               lift (twist f xs, unit f `cons` fs))
    $ initial


-- | Transformation of arbitrary length base on Bluestein on a 5-smooth size.
--
transformChirp235
    :: (Shape sh, Slice sh, Numeric e)
    => Mode
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
transformChirp235 mode arr =
  let n = indexHead (shape arr)
      f = ceiling5Smooth (2*n)
  in
  transformChirp mode f (dit235 Forward) (dit235 Inverse) arr


transformChirp
    :: (Shape sh, Slice sh, Numeric e)
    => Mode
    -> Exp Int
    -> (forall sh'. (Shape sh', Slice sh') => Acc (Array (sh':.Int) (Complex e)) -> Acc (Array (sh':.Int) (Complex e)))
    -> (forall sh'. (Shape sh', Slice sh') => Acc (Array (sh':.Int) (Complex e)) -> Acc (Array (sh':.Int) (Complex e)))
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
transformChirp mode p analysis synthesis arr =
  let sz :. n   = unlift (shape arr)
      --
      chirp     =
        generate (index1 p) $ \ix ->
          let k  = unindex1 ix
              sk = fromIntegral (if p > 2*k then k else k-p)
              w  = pi * sk * sk / fromIntegral n
          in
          lift $ cos w :+ signOfMode mode * sin w
      --
      spectrum  = analysis
                $ map conjugate chirp
                  `consV`
                  reshape (lift (Z :. shapeSize sz :. p))
                          (pad p 0 (zipWithExtrude1 (*) chirp arr))
      scaleDown xs =
        let scale x (unlift -> r :+ i) = lift (x*r :+ x*i)
            len                        = indexHead (shape xs)
        in  map (scale (recip (fromIntegral len))) xs
  in
  if n <= 1
    then arr
    else take n
       $ scaleDown
       $ zipWithExtrude1 (*) chirp
       $ synthesis
       $ zipWithExtrude1 (*) (headV spectrum)
       $ reshape (lift (sz:.p)) (tailV spectrum)


ceiling5Smooth :: Exp Int -> Exp Int
ceiling5Smooth n =
  let (i2,i3,i5) = unlift (snd (ceiling5Smooth' (fromIntegral n :: Exp Double)))
  in  pow i2 2 * pow i3 3 * pow i5 5

ceiling5Smooth'
    :: (RealFloat a, Ord a, FromIntegral Int a)
    => Exp a
    -> Exp (a, (Int,Int,Int))
ceiling5Smooth' n =
  let d3 = ceiling (logBase 3 n)
      d5 = ceiling (logBase 5 n)
      --
      argmin x y = if fst x < fst y then x else y
  in
  the $ fold1All argmin
      $ generate (index2 d5 d3) -- this is probably quite small!
                 (\(unlift -> Z :. i5 :. i3) ->
                    let
                        p53 = 5 ** fromIntegral i5 * 3 ** fromIntegral i3
                        i2  = 0 `max` ceiling (logBase 2 (n/p53))
                    in
                    lift ( p53 * 2 ** fromIntegral i2
                         , (i2,i3,i5)
                         ))

-- Utilities
-- ---------

pow :: Exp Int -> Exp Int -> Exp Int
pow x k
  = snd
  $ while (\ip -> fst ip < k)
          (\ip -> lift (fst ip + 1, snd ip * x))
          (lift (0::Int,1::Int))

pad :: (Shape sh, Slice sh, Elt e)
    => Exp Int
    -> Exp e
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
pad n x xs =
  let sz = indexTail (shape xs)
      sh = lift (sz :. n)
  in
  xs ++ fill sh x

cons :: forall sh e. (Shape sh, Slice sh, Elt e)
     => Acc (Array sh e)
     -> Acc (Array (sh:.Int) e)
     -> Acc (Array (sh:.Int) e)
cons x xs =
  let x' = reshape (lift (shape x :. (1::Int))) x
  in  x' ++ xs

consV :: forall sh e. (Shape sh, Slice sh, Elt e)
      => Acc (Array (sh:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
consV x xs =
  let sh :. n = unlift (shape x) :: Exp sh :. Exp Int
      newsh = lift (sh :. (1::Int) :. n) :: Exp (sh :. Int :. Int)
  in  reshape newsh x ++^ xs

-- If input is regular, stays regular
headV :: (Shape sh, Slice sh, Elt e)
      => Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int) e)
headV xs = slice xs (lift (Any :. (0 :: Exp Int) :. All))

tailV :: forall sh e. (Shape sh, Slice sh, Elt e)
      => Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
tailV = tailOn _2

-- If Exp Int and array shape are same, then stays regular.
dropV :: forall sh e. (Shape sh, Slice sh, Elt e)
      => Exp Int
      -> Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
dropV = dropOn _2

--New shape is fully determined by old one and fac. If shape and fac same, stays regular.
sieve
    :: forall sh e. (Shape sh, Slice sh, Elt e)
    => Exp Int
    -> Exp Int
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
sieve fac start xs =
  let sh :. n = unlift (shape xs) :: Exp sh :. Exp Int
  in
  backpermute
    (lift (sh :. n `quot` fac))
    (\(unlift -> ix :. j :: Exp sh :. Exp Int) -> lift (ix :. fac*j + start))
    xs

--New shape is fully determined by old one and fac. If shape and fac same, stays regular.
sieveV
    :: forall sh e. (Shape sh, Slice sh, Elt e)
    => Exp Int
    -> Exp Int
    -> Acc (Array (sh:.Int:.Int) e)
    -> Acc (Array (sh:.Int:.Int) e)
sieveV fac start xs =
  let sh :. m :. n = unlift (shape xs) :: Exp sh :. Exp Int :. Exp Int
  in
  backpermute
    (lift (sh :. m `quot` fac :. n))
    (\(unlift -> ix :. j :. i :: Exp sh :. Exp Int :. Exp Int) -> lift (ix :. fac*j+start :. i))
    xs

-- If fac and shape stay the same, then twist stays regular.
twist :: forall sh e. (Shape sh, Slice sh, Elt e)
      => Exp Int
      -> Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
twist fac xs =
  let sh :. m :. n = unlift (shape xs) :: Exp sh :. Exp Int :. Exp Int
  in
  backpermute
    (lift (sh :. fac*m :. n `quot` fac))
    (\(unlift -> ix :. j :. i :: Exp sh :. Exp Int :. Exp Int) -> lift (ix :. j `quot` fac :. fac*i + j `rem` fac))
    xs


infixr 5 ++^
(++^) :: forall sh e. (Shape sh, Slice sh, Elt e)
      => Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
(++^) = concatOn _2

-- if xs and ys regulary shaped, then zipWithExtrude1 also regulary shaped
zipWithExtrude1
    :: (Shape sh, Slice sh, Elt a, Elt b, Elt c)
    => (Exp a -> Exp b -> Exp c)
    -> Acc (Array DIM1      a)
    -> Acc (Array (sh:.Int) b)
    -> Acc (Array (sh:.Int) c)
zipWithExtrude1 f xs ys =
  zipWith f (replicate (lift (indexTail (shape ys) :. All)) xs) ys

zipWithExtrude2
    :: (Shape sh, Slice sh, Elt a, Elt b, Elt c)
    => (Exp a -> Exp b -> Exp c)
    -> Acc (Array DIM2           a)
    -> Acc (Array (sh:.Int:.Int) b)
    -> Acc (Array (sh:.Int:.Int) c)
zipWithExtrude2 f xs ys =
  zipWith f (replicate (lift (indexTail (indexTail (shape ys)) :. All :. All)) xs) ys

transpose
    :: forall sh e. (Shape sh, Slice sh, Elt e)
    => Acc (Array (sh:.Int:.Int) e)
    -> Acc (Array (sh:.Int:.Int) e)
transpose = transposeOn _1 _2

-- Transform2 stays regular, fully based on shape (and constant expresion)
transform2
    :: (Shape sh, Slice sh, Num e)
    => Exp e
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
transform2 v xs =
  generate
    (lift (indexTail (shape xs) :. (2::Int)))
    (\(unlift -> ix :. k :: Exp sh :. Exp Int) ->
        let x0 = xs ! lift (ix :. (0::Int))
            x1 = xs ! lift (ix :. (1::Int))
        in
        if k == 0 then x0+x1
                  else x0+v*x1)

transform3
    :: forall sh e. (Shape sh, Slice sh, Num e)
    => Exp (e,e)
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
transform3 (unlift -> (z1,z2)) xs =
  generate
    (lift (indexTail (shape xs) :. (3::Int)))
    (\(unlift -> ix :. k :: Exp sh :. Exp Int) ->
        let
            x0 = xs ! lift (ix :. (0::Int))
            x1 = xs ! lift (ix :. (1::Int))
            x2 = xs ! lift (ix :. (2::Int))
            --
            ((s,_), (zx1,zx2)) = sumAndConvolve2 (x1,x2) (z1,z2)
        in
        if k == 0    then x0 + s   else
        if k == 1    then x0 + zx1
        {- k == 2 -} else x0 + zx2)

transform4
    :: forall sh e. (Shape sh, Slice sh, Num e)
    => Exp (e,e,e)
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
transform4 (unlift -> (z1,z2,z3)) xs =
  generate
    (lift (indexTail (shape xs) :. (4::Int)))
    (\(unlift -> ix :. k :: Exp sh :. Exp Int) ->
        let
            x0 = xs ! lift (ix :. (0::Int))
            x1 = xs ! lift (ix :. (1::Int))
            x2 = xs ! lift (ix :. (2::Int))
            x3 = xs ! lift (ix :. (3::Int))
            --
            x02a = x0+x2
            x02b = x0+z2*x2
            x13a = x1+x3
            x13b = x1+z2*x3
        in
        if k == 0    then x02a +      x13a else
        if k == 1    then x02b + z1 * x13b else
        if k == 2    then x02a + z2 * x13a
        {- k == 3 -} else x02b + z3 * x13b)

-- Use Rader's trick for mapping the transform to a convolution and apply
-- Karatsuba's trick at two levels (i.e. total three times) to that convolution.
--
-- 0 0 0 0 0
-- 0 1 2 3 4
-- 0 2 4 1 3
-- 0 3 1 4 2
-- 0 4 3 2 1
--
-- Permutation.T: 0 1 2 4 3
--
-- 0 0 0 0 0
-- 0 1 2 4 3
-- 0 2 4 3 1
-- 0 4 3 1 2
-- 0 3 1 2 4
--
transform5
    :: forall sh e. (Shape sh, Slice sh, Num e)
    => Exp (e,e,e,e)
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
transform5 (unlift -> (z1,z2,z3,z4)) xs =
  generate
    (lift (indexTail (shape xs) :. (5::Int)))
    (\(unlift -> ix :. k :: Exp sh :. Exp Int) ->
        let
            x0 = xs ! lift (ix :. (0::Int))
            x1 = xs ! lift (ix :. (1::Int))
            x2 = xs ! lift (ix :. (2::Int))
            x3 = xs ! lift (ix :. (3::Int))
            x4 = xs ! lift (ix :. (4::Int))
            --
            ((s,_), (d1,d2,d4,d3)) = sumAndConvolve4 (x1,x3,x4,x2) (z1,z2,z4,z3)
        in
        if k == 0    then x0 + s  else
        if k == 1    then x0 + d1 else
        if k == 2    then x0 + d2 else
        if k == 3    then x0 + d3
        {- k == 4 -} else x0 + d4)


-- Some small size convolutions using the Karatsuba trick.
--
-- This does not use Toom-3 multiplication, because this requires division by
-- 2 and 6, and thus 'Fractional' constraints.
--
sumAndConvolve2
    :: Num e
    => (Exp e, Exp e)
    -> (Exp e, Exp e)
    -> ((Exp e, Exp e), (Exp e, Exp e))
sumAndConvolve2 (a0,a1) (b0,b1) =
  let sa01   = a0+a1
      sb01   = b0+b1
      ab0ab1 = a0*b0+a1*b1
  in
  ((sa01, sb01), (ab0ab1, sa01*sb01-ab0ab1))

-- sumAndConvolve3
--     :: Num e
--     => (Exp e, Exp e, Exp e)
--     -> (Exp e, Exp e, Exp e)
--     -> ((Exp e, Exp e), (Exp e, Exp e, Exp e))
-- sumAndConvolve3 (a0,a1,a2) (b0,b1,b2) =
--   let ab0   = a0*b0
--       dab12 = a1*b1 - a2*b2
--       sa01  = a0+a1; sb01 = b0+b1; tab01 = sa01*sb01 - ab0
--       sa02  = a0+a2; sb02 = b0+b2; tab02 = sa02*sb02 - ab0
--       sa012 = sa01+a2
--       sb012 = sb01+b2
--       --
--       d0    = sa012*sb012 - tab01 - tab02
--       d1    = tab01 - dab12
--       d2    = tab02 + dab12
--   in
--   ((sa012, sb012), (d0, d1, d2))

sumAndConvolve4
  :: Num e
  => (Exp e, Exp e, Exp e, Exp e)
  -> (Exp e, Exp e, Exp e, Exp e)
  -> ((Exp e, Exp e), (Exp e, Exp e, Exp e, Exp e))
sumAndConvolve4 (a0,a1,a2,a3) (b0,b1,b2,b3) =
  let ab0    = a0*b0
      ab1    = a1*b1
      sa01   = a0+a1; sb01 = b0+b1
      ab01   = sa01*sb01 - (ab0+ab1)
      ab2    = a2*b2
      ab3    = a3*b3
      sa23   = a2+a3; sb23 = b2+b3
      ab23   = sa23*sb23 - (ab2+ab3)
      c0     = ab0  + ab2 - (ab1 + ab3)
      c1     = ab01 + ab23
      ab02   = (a0+a2)*(b0+b2)
      ab13   = (a1+a3)*(b1+b3)
      sa0123 = sa01+sa23
      sb0123 = sb01+sb23
      ab0123 = sa0123*sb0123 - (ab02+ab13)
      --
      d0     = ab13   + c0
      d1     = c1
      d2     = ab02   - c0
      d3     = ab0123 - c1
  in
  ((sa0123, sb0123), (d0, d1, d2, d3))

makeExp :: Numeric e => Int -> Exp e
makeExp = fromIntegral . constant

signOfMode :: P.Num a => Mode -> a
signOfMode m
  = case m of
      Forward -> -1
      Reverse ->  1
      Inverse ->  1





------------------------------------------
-- For reference the older version, maybe compiling this is faster..

-- Rank-generalised Cooley-Tuckey DFT
--
-- We require the innermost dimension be passed as a Haskell value because we
-- can't do divide-and-conquer recursion directly in the meta-language.
--
fftV2 :: forall sh e. (Slice sh, Shape sh, A.RealFloat e, A.FromIntegral Int e)
    => e
    -> sh
    -> Int
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
fftV2 sign sh sz arr
  | P.any (P.not . isPow2) (shapeToList (sh:.sz))
  = error $ printf "fft: array dimensions must be powers-of-two, but are: %s" (showShape (sh:.sz))
  --
  | P.otherwise
  = go sz 0 1
  where
    go :: Int -> Int -> Int -> Acc (Array (sh:.Int) (Complex e))
    go len offset stride
      | len P.== 2
      = A.generate (constant (sh :. len)) swivel

      | P.otherwise
      = combine
          (go (len `div` 2) offset            (stride * 2))
          (go (len `div` 2) (offset + stride) (stride * 2))

      where
        len'    = the (unit (constant len))
        offset' = the (unit (constant offset))
        stride' = the (unit (constant stride))

        swivel ix =
          let sh' :. sz' = unlift ix :: Exp sh :. Exp Int
          in
          sz' A.== 0 ? ( (arr ! lift (sh' :. offset')) + (arr ! lift (sh' :. offset' + stride'))
          {-  A.== 1-} , (arr ! lift (sh' :. offset')) - (arr ! lift (sh' :. offset' + stride')) )

        combine evens odds =
          let odds' = A.generate (A.shape odds) (\ix -> twiddle len' (indexHead ix) * odds!ix)
          in
          append (A.zipWith (+) evens odds') (A.zipWith (-) evens odds')

        twiddle n' i' =
          let n = A.fromIntegral n'
              i = A.fromIntegral i'
              k = 2*pi*i/n
          in
          lift ( cos k :+ A.constant sign * sin k )


-- Append two arrays. This is a specialised version of (A.++) which does not do
-- bounds checking or intersection.
--
append
    :: forall sh e. (Slice sh, Shape sh, Elt e)
    => Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
append xs ys
  = let sh :. n = unlift (A.shape xs)     :: Exp sh :. Exp Int
        _  :. m = unlift (A.shape ys)     :: Exp sh :. Exp Int
    in
    generate (lift (sh :. n+m))
             (\ix -> let sz :. i = unlift ix :: Exp sh :. Exp Int
                     in  i A.< n ? (xs ! lift (sz:.i), ys ! lift (sz:.i-n) ))


isPow2 :: Int -> Bool
-- isPow2 0 = True
-- isPow2 1 = False
isPow2 x | x P.== 0    = True
         | x P.== 1    = False
         | P.otherwise = x B..&. (x P.- 1) P.== 0
