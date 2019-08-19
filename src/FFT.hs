{-# language FlexibleContexts    #-}
{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE PatternGuards       #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
{-# LANGUAGE TypeFamilies        #-}
{-# LANGUAGE TypeOperators       #-}

module FFT where

import Types
import Utility

import Data.Array.Accelerate                                        as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                              as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                           as A
import Data.Array.Accelerate.Math.DFT.Centre                        as A
import Data.Array.Accelerate.Math.FFT                               as A 
import Data.Array.Accelerate.Math.FFT.Type                          as A
--import Data.Array.Accelerate.Math.FFT.LLVM.Native.Ix
--import Data.Array.Accelerate.Math.FFT.LLVM.Native.Base

import Data.Array.Accelerate.LLVM.Native                            as CPU
import Data.Array.Accelerate.Array.Sugar                            as S hiding (shape)
import Data.Array.Accelerate.LLVM.Native.Foreign
import Data.Array.Accelerate.Array.Unique                           as A
import Data.Array.Accelerate.Array.Lifted                           as A
import Data.Array.Accelerate.Error                                  as A
import Data.Array.Accelerate.Type                                   as A
import Data.Array.Accelerate.Lifetime                               as A


import Data.Ix                                                      ( Ix )
--import Data.Array.CArray                                            ( CArray(..) )
import Data.Array.CArray.Base                                       ( CArray(..) )
import qualified Data.Array.CArray                                  as C
import qualified Data.Array.CArray.Base                             as C
import Foreign.Ptr
import Foreign.Storable
import Foreign.ForeignPtr

import Math.FFT                                                     as FFT
import Math.FFT.Base                                                ( FFTWReal, Sign(..), Flag, measure, destroyInput )
import Text.Printf
import Data.Bits

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..), maybe, fail, (=<<))
import Debug.Trace

----------------------
-- Fourier transformations
fft :: Acc (Matrix Visibility) -> Acc (Matrix Visibility)
fft = shift2D . myfft2D Forward . ishift2D

fftExpand :: Acc (Matrix Visibility) -> Acc (Matrix Visibility)
fftExpand m = extract_mid n_ . shift2D . myfft2D Forward . ishift2D . (`pad_mid` n) $ m
    where
        (Z :. n_ :. _) = (unlift . shape) m :: Z :. Exp Int :. Exp Int
        pw = A.ceiling (logBase 2 (A.fromIntegral n_) :: Exp F) :: Exp Int
        n = 2 ^ pw :: Exp Int


ifft :: Acc (Matrix Visibility) -> Acc (Matrix Visibility)
ifft = shift2D . myfft2D Inverse . ishift2D

ifftExpand :: Acc (Matrix Visibility) -> Acc (Matrix Visibility)
ifftExpand m = extract_mid n_ . shift2D . myfft2D Inverse . ishift2D . (`pad_mid` n) $ m
    where
        (Z :. n_ :. _) = (unlift . shape) m :: Z :. Exp Int :. Exp Int
        pw = A.ceiling (logBase 2 (A.fromIntegral n_) :: Exp F) :: Exp Int
        n = 2 ^ pw :: Exp Int


myfft2D :: Numeric e => Mode -> Acc (Array DIM2 (Complex e)) -> Acc (Array DIM2 (Complex e))
myfft2D = fft2D --myfft2DFor
      
myfft2DFor :: Numeric e => Mode -> Acc (Array DIM2 (Complex e)) -> Acc (Array DIM2 (Complex e))
myfft2DFor mode = foreignAcc (fft2DVect mode) $ A.map (\_ -> -89978978.4e0100) {-Bollocks implementation, for checking-}
    where
        fft2DAvoid :: forall e. (Numeric e) 
                    => Mode 
                    -> ForeignAcc (Array DIM2 (Complex e) -> Array DIM2 (Complex e))
        fft2DAvoid mode = ForeignAcc (nameOf mode (undefined::DIM2))
            $ case numericR::NumericR e of
              NumericRfloat32 -> liftIO . liftAtoC go
              NumericRfloat64 -> liftIO . liftAtoC go
            where
                go :: FFTWReal r => CArray (Int,Int) (Complex r) -> CArray (Int,Int) (Complex r)
                go = FFT.dftG (signOf mode) flags [0,1]

        fft2DRegular :: forall e. (Numeric e) 
                    => Mode 
                    -> ForeignAcc (Array DIM3 (Complex e) -> Array DIM3 (Complex e))
        fft2DRegular mode = ForeignAcc (nameOfV mode (undefined::DIM2))
            $ case numericR::NumericR e of
              NumericRfloat32 -> liftIO . liftAtoC go
              NumericRfloat64 -> liftIO . liftAtoC go
            where
                go :: (P.Show r, FFTWReal r) => CArray (Int,Int,Int) (Complex r) -> CArray (Int,Int,Int) (Complex r)
                go = FFT.dftG (signOf mode) flags [1,2]

        fft2DVect :: forall e. (Numeric e)
                => Mode -> VectorisedForeign (Array DIM2 (Complex e) -> Array DIM2 (Complex e))
        fft2DVect mode = VectorisedForeign $  f
            where
                f :: Arrays a' => LiftedType (Array DIM2 (Complex e)) a' -> LiftedType (Array DIM2 (Complex e)) b' -> ForeignAcc (a' -> b')
                f AvoidedT AvoidedT = fft2DAvoid mode
                f RegularT RegularT = fft2DRegular mode
                f IrregularT IrregularT = error "no irregular stuff"

myfft2D32 :: Numeric e => Mode -> Acc (Array DIM2 (Complex e)) -> Acc (Array DIM2 (Complex e))
myfft2D32 = fft2D

signOf :: Mode -> Sign
signOf Forward = DFTForward
signOf _       = DFTBackward

flags :: Flag
flags = estimate .|. destroyInput

nameOf :: forall sh. Shape sh => Mode -> sh -> String
nameOf Forward _ = printf "FFTW.dft%dD"  (rank (undefined::sh))
nameOf _       _ = printf "FFTW.idft%dD" (rank (undefined::sh))

nameOfV :: forall sh. Shape sh => Mode -> sh -> String
nameOfV Forward _ = printf "FFTW.dft%dDV"  (rank (undefined::sh))
nameOfV _       _ = printf "FFTW.idft%dDV" (rank (undefined::sh))

-- | Lift an operation on CArray into an operation on Accelerate arrays
--
liftAtoC
    :: (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Shape sh, Ix ix, Elt ix, Elt e, IsFloating e, Storable e', ArrayPtrs e ~ Ptr e')
    => (CArray ix (Complex e') -> CArray ix (Complex e'))
    -> Array sh (Complex e)
    -> IO (Array sh (Complex e))
liftAtoC f a = c2a . f =<< a2c a


-- | Convert a multidimensional Accelerate array of complex numbers into
-- a packed CArray of complex numbers suitable for use by FFTW.
--
a2c :: forall ix sh e e'. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Ix ix, Elt ix, Shape sh, IsFloating e, Storable e', ArrayPtrs e ~ Ptr e')
    => Array sh (Complex e)
    -> IO (CArray ix (Complex e'))
a2c arr
  | FloatingDict <- floatingDict (floatingType :: FloatingType e)
  = let
        (lo,hi) = shapeToRange (arrayShape arr)
        bnds    = (fromIxShapeRepr lo, fromIxShapeRepr hi)
        n       = S.size (arrayShape arr)
    in
    C.createCArray       bnds $ \p_cs      ->
    withComplexArrayPtrs arr  $ \p_re p_im ->
      let
          -- TLM: Should we execute this in parallel using the worker threads of
          -- the current target? (Native)
          go !i | i P.>= n = return ()
          go !i            = do
            re <- peekElemOff p_re i
            im <- peekElemOff p_im i
            pokeElemOff p_cs i (re :+ im)
            go (i+1)
      in
      go 0


-- | Convert a packed CArray of complex numbers into an unzipped (SoA)
-- multidimensional Accelerate array of complex numbers.
--
c2a :: forall ix sh e e'. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Ix ix, Elt ix, Shape sh, Elt e, IsFloating e, Storable e', ArrayPtrs e ~ Ptr e')
    => CArray ix (Complex e')
    -> IO (Array sh (Complex e))
c2a carr
  | FloatingDict <- floatingDict (floatingType :: FloatingType e)
  = let
        (lo,hi) = C.bounds carr
        n       = C.rangeSize (lo,hi)
        sh      = rangeToShape (toIxShapeRepr lo, toIxShapeRepr hi)
    in do
      arr <- allocateArray sh
      C.withCArray carr        $ \p_cs      ->
        withComplexArrayPtrs arr $ \p_re p_im -> do
          let go !i | i P.>= n = return ()
              go !i            = do
                  re :+ im <- peekElemOff p_cs i
                  pokeElemOff p_re i re
                  pokeElemOff p_im i im
                  go (i+1);
          --
          go 0;
          return arr;



-- Dig out the underlying pointers of the Accelerate SoA data structure
--

withComplexArrayPtrs
    :: forall sh e a. IsFloating e
    => Array sh (Complex e)
    -> (ArrayPtrs e -> ArrayPtrs e -> IO a)
    -> IO a
withComplexArrayPtrs (Array _ adata) k
  | AD_Pair (AD_Pair AD_Unit ad1) ad2 <- adata
  = case floatingType :: FloatingType e of
      TypeFloat{}   -> withArrayData arrayElt ad1 $ \p1 -> withArrayData arrayElt ad2 $ \p2 -> k p1 p2
      TypeDouble{}  -> withArrayData arrayElt ad1 $ \p1 -> withArrayData arrayElt ad2 $ \p2 -> k p1 p2
      TypeCFloat{}  -> withArrayData arrayElt ad1 $ \p1 -> withArrayData arrayElt ad2 $ \p2 -> k p1 p2
      TypeCDouble{} -> withArrayData arrayElt ad1 $ \p1 -> withArrayData arrayElt ad2 $ \p2 -> k p1 p2

-- withScalarArrayPtrs
--     :: forall sh e a. IsFloating e
--     => Array sh e
--     -> (ArrayPtrs e -> IO a)
--     -> IO a
-- withScalarArrayPtrs (Array _ adata) =
--   case floatingType :: FloatingType e of
--     TypeFloat{}   -> withArrayData arrayElt adata
--     TypeDouble{}  -> withArrayData arrayElt adata
--     TypeCFloat{}  -> withArrayData arrayElt adata
--     TypeCDouble{} -> withArrayData arrayElt adata

withArrayData
    :: (ArrayPtrs e ~ Ptr a)
    => ArrayEltR e
    -> ArrayData e
    -> (Ptr a -> IO b)
    -> IO b
withArrayData ArrayEltRfloat   (AD_Float   ua) = withUniqueArrayPtr ua
withArrayData ArrayEltRdouble  (AD_Double  ua) = withUniqueArrayPtr ua
withArrayData ArrayEltRcfloat  (AD_CFloat  ua) = withUniqueArrayPtr ua
withArrayData ArrayEltRcdouble (AD_CDouble ua) = withUniqueArrayPtr ua
withArrayData _ _ =
  $internalError "withArrayData" "expected array of [C]Float or [C]Double"


-- Converting between Accelerate multidimensional shapes/indices and those used
-- by the CArray package (Data.Ix)
--

type family IxShapeRepr e where
  IxShapeRepr ()    = ()
  IxShapeRepr Int   = ((),Int)
  IxShapeRepr (t,h) = (IxShapeRepr t, h)

fromIxShapeRepr
    :: forall ix sh. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Shape sh, Elt ix)
    => sh
    -> ix
fromIxShapeRepr = liftToElt (go (eltType (undefined::ix)))
  where
    go :: forall ix'. TupleType ix' -> IxShapeRepr ix' -> ix'
    go UnitTuple                                                 ()     = ()
    go (PairTuple tt _)                                          (t, h) = (go tt t, h)
    go (SingleTuple (NumScalarType (IntegralNumType TypeInt{}))) ((),h) = h
    go _ _
      = $internalError "fromIxShapeRepr" "expected Int dimensions"

toIxShapeRepr
    :: forall ix sh. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Shape sh, Elt ix)
    => ix
    -> sh
toIxShapeRepr = liftToElt (go (eltType (undefined::ix)))
  where
    go :: forall ix'. TupleType ix' -> ix' -> IxShapeRepr ix'
    go UnitTuple        ()                                             = ()
    go (SingleTuple     (NumScalarType (IntegralNumType TypeInt{}))) h = ((), h)
    go (PairTuple tt _) (t, h)                                         = (go tt t, h)
    go _ _
      = error "toIxShapeRepr: not a valid Data.Ix index"


{-
{-# INLINE liftCtoA #-}
liftCtoA
    :: forall ix sh e. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Shape sh, Elt ix, Numeric e)
    => (CArray ix (Complex e) -> CArray ix (Complex e))
    -> Array sh (Complex e)
    -> IO (Array sh (Complex e))
liftCtoA f a =
  newFull =<< liftIO (withCArray a (fromCArray . f))

-- /O(1)/ Convert a CArray to an Accelerate array
--
{-# INLINE fromCArray #-}
fromCArray
    :: forall ix sh e. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Shape sh, Elt ix, Numeric e)
    => CArray ix (Complex e)
    -> IO (Array sh (Complex e))
fromCArray (CArray lo hi _ fp) = do
  --
  sh <- return $ rangeToShape (toIxShapeRepr lo, toIxShapeRepr hi) :: IO sh
  ua <- newUniqueArray (castForeignPtr fp :: ForeignPtr e)
  --
  case numericR::NumericR e of
    NumericRfloat32 -> return $ Array (fromElt sh) (AD_Vec 2# (AD_Float  ua))
    NumericRfloat64 -> return $ Array (fromElt sh) (AD_Vec 2# (AD_Double ua))

-- /O(1)/ Use an Accelerate array as a CArray
--
{-# INLINE withCArray #-}
withCArray
    :: forall ix sh e a. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Shape sh, Elt ix, Numeric e)
    => Array sh (Complex e)
    -> (CArray ix (Complex e) -> IO a)
    -> IO a
withCArray arr f =
  let
      sh        = shape arr
      (lo, hi)  = shapeToRange sh
      wrap fp   = CArray (fromIxShapeRepr lo) (fromIxShapeRepr hi) (S.size sh) (castForeignPtr fp)
  in
  withArray arr (f . wrap)


-- Use underlying array pointers
--
{-# INLINE withArray #-}
withArray
    :: forall sh e a. Numeric e
    => Array sh (Complex e)
    -> (ForeignPtr e -> IO a)
    -> IO a
withArray (Array _ adata) = withArrayData (numericR::NumericR e) adata

{-# INLINE withArrayData #-}
withArrayData
    :: NumericR e
    -> ArrayData (EltRepr (Complex e))
    -> (ForeignPtr e -> IO a)
    -> IO a
withArrayData NumericRfloat32 ((AD_Float  ua)) = withLifetime (uniqueArrayData ua)
withArrayData NumericRfloat64 ((AD_Double ua)) = withLifetime (uniqueArrayData ua)
--withArrayData NumericRfloat32 (AD_Vec _ (AD_Float  ua)) = withLifetime (uniqueArrayData ua)
--withArrayData NumericRfloat64 (AD_Vec _ (AD_Double ua)) = withLifetime (uniqueArrayData ua)
-}