module PermuteSeqError where
import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)

import Data.Array.Accelerate.LLVM.Native                  as CPU

testF :: Acc (Scalar Int) -> Acc (Array DIM2 Int)
testF i' = let
    i = the i'
    orig = generate (index2 4 4) shapeSize
    def = fill (index2 5 5) i
 in permute const def id orig

permuteTest3 = collect . elements $ produce 3 testF