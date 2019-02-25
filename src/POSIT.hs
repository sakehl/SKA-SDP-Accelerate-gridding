{-# LANGUAGE DeriveGeneric         #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE PatternSynonyms       #-}
{-# LANGUAGE RebindableSyntax      #-}
{-# LANGUAGE ScopedTypeVariables   #-}
{-# LANGUAGE ViewPatterns          #-}

module POSIT (

  Posit16(..), pattern P16_,
  p16_to_f32, f32_to_p16,

) where

import Data.Array.Accelerate                    as A
import Data.Array.Accelerate.Data.Bits          as A hiding ( bit, shift )


newtype Posit16 = P16 Word16
  deriving (Show, Generic)

instance Elt Posit16
instance IsProduct Elt Posit16

pattern P16_ :: Exp Word16 -> Exp Posit16
pattern P16_ p = Pattern p
{-# COMPLETE P16_ #-}


infixl 8 .<<.
(.<<.) :: Bits a => Exp a -> Exp Int -> Exp a
(.<<.) = unsafeShiftL

infixl 8 .>>.
(.>>.) :: Bits a => Exp a -> Exp Int -> Exp a
(.>>.) = unsafeShiftR

-- Convert 16-bit Posit to 32-bit IEEE Float
--
p16_to_f32 :: Exp Posit16 -> Exp Float
p16_to_f32 (P16_ p) = bitcast f
  where
    sign :: Exp Word32
    sign = if p .&. 0x8000 /= 0 then 0x80000000 else 0
    f    = if p .&. 0x7FFF /= 0
             -- decode |p|
             then
               let -- decode regime
                   f0 = sign /= 0 ? ( 0x10000 - p, p )

                   u  = f0 .&. 0x4000 /= 0
                   v  = u ? (complement f0, f0)
                   w  = countLeadingZeros (v .<<. 2)
                   x  = 0x1000000 * w

                   f1 = fromIntegral f0 .<<. w  -- XXX: requires 32-bits
                   s1 = u ? ( 0x3f800000 + x    -- XXX: must be signed
                            , 0x3e800000 - x)

                   -- decode exponent bit
                   s2 =
                     if f1 .&. 0x1000 /= 0
                       then s1 + 0x800000
                       else s1
               in
               sign .|. fromIntegral s2 .|. ((f1 .&. 0x0FFF) .<<. 11)

             -- exception cases NaN and zero
             else
               sign /= 0 ? ( 0x7fffffff, 0)


-- Convert 32-bit IEEE Float to 16-bit Posit with correct rounding
--
f32_to_p16 :: Exp Float -> Exp Posit16
f32_to_p16 (bitcast -> f :: Exp Word32) = P16_ (fromIntegral p)
  where
    p = if f .&. 0x7f800000 == 0x7f800000
          -- ±infinity, NaN become NaR
          then 0x8000
          else
            -- work with |f|
            let p0 = f .&. 0x7fffffff in
            if  p0 == 0
              -- ±0 both become 0
              then 0
              else
                let p1 = if p0 > 0x4d000000 then 0x7fff else    -- |f| > 2^27 rounds to maxpos
                         if p0 < 0x32000000 then 1      else    -- |f| < 2^-27 rounds to minpos
                           let
                             s0 = (fromIntegral p0 .>>. 23) - 127 -- power-of-2 meaning of float exponent field
                             b0 = s0 .&. 1                        -- posit exponent bit

                             -- determine number of extra regime bits needed
                             T2 s2 p2 =
                               if s0 >= 0
                                 then let s1 = s0 .>>. 1
                                      in  T2 s1 (0x3ffffff - (0x1ffffff .>>. s1))
                                 else let s1 = (-s0-1) .>>. 1
                                      in  T2 s1 (0x1000000 .>>. s1)

                             -- install exponent bit
                             p3 = if b0 /= 0
                                    then p2 .|. (0x800000 .>>. s2)
                                    else p2

                             -- unrounded fraction
                             p4 :: Exp Word64       -- XXX: requires 64-bits
                             p4 = (p3 .<<. s2) .|. (fromIntegral f .&. 0x7fffff)

                             b1 = 0x400 .<<. s2 -- bit n+1
                             p5 = if p4 .&. b1 /= 0 && (p4 .&. (b1-1) /= 0 || p4 .&. (b1 .<<. 1) /= 0)
                                    then p4 + b1
                                    else p4
                           in
                           p5 .>>. (s2 + 11)
                in
                testBit f 31 ? ( 0x10000 - p1, p1 )

