module Types where

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A

--Floating point
type F = Double
--Integer
type I = Int

type BaseLine = F
type BaseLines = (F,F,F)
type Antenna = Int
type Visibility = Complex F
type Time = F
type Frequency = F
type Image = Matrix F
type FourierSpace = Matrix Visibility

type Kernel = Array DIM4 Visibility

type WKernel = Array DIM4 Visibility
type AKernel = Array DIM2 Visibility
type AWKernel = Array DIM2 Visibility