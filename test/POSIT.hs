{-# LANGUAGE RebindableSyntax #-}
{-# LANGUAGE TypeApplications #-}

import POSIT

import Data.Array.Accelerate             as A
import Data.Array.Accelerate.Debug       as A
import Data.Array.Accelerate.Data.Bits   as A
import Data.Array.Accelerate.Interpreter as I
import Data.Array.Accelerate.LLVM.Native as CPU


roundtrip :: Exp Word16 -> Exp Word16
roundtrip p =
  let f      = p16_to_f32 (P16_ p)
      P16_ q = f32_to_p16 f
  in
  q

exhaustive :: Acc (Vector (Word16, Word16), Scalar Int)
exhaustive =
  let ws  = enumFromN (index1 65536) 0
      ws' = map roundtrip ws
  in
  filter (uncurry (/=)) (zip ws ws')

