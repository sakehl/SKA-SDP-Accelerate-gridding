{-# language FlexibleContexts    #-}
{-# language RebindableSyntax    #-}
{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}

module Main where

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A
import Data.Array.Accelerate.Math.DFT.Centre              as A hiding (shift1D, shift2D, shift3D)
import Data.Array.Accelerate.Math.FFT                     as A

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I

import Gridding

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..))
import Data.Char (isSpace)
import Text.Printf
import Data.List (intercalate)

import Debug.Trace
import System.IO.Unsafe (unsafePerformIO)
import Control.Exception (assert)

import Control.Lens as L (_1, _2)

main :: IO ()
main = testWCache

----------------------
-- I/O things (for testing)
testSimple :: IO ()
testSimple = do
    testData <- testData0
    let n   = (P.length . P.fst) testData
        m   = (P.length . P.snd) testData
        size = assert (n P.== m) $ Z :. n
        vis = fromList size $ P.fst testData
        uvw = fromList size $ P.snd testData
        src      = undefined
        theta    = 2*0.05
        lam      = 18000
        kwargs   = undefined
        otargs   = undefined
        (d,p, _) = CPU.run $ do_imaging theta lam uvw src vis simple_imaging kwargs otargs
        (muvw, mvis) = unzip $ mirror_uvw (use uvw) (use vis)

        fst i j = indexArray d (Z :. i :. j)
        maxi = P.maximum (toList d)
    P.writeFile "data/mirror_uvw.csv" (makeVFile (showTriple $ printf "%e") (CPU.run muvw))
    P.writeFile "data/result.csv" (makeMFile (printf "%e") d)
    P.writeFile "data/result_p.csv" (makeMFile (printf "%e") p)
    -- P.putStrLn (P.show maxi)

testWCache :: IO ()
testWCache = do
    testData <- testData0
    let n   = (P.length . P.fst) testData
        n2   = (P.length . P.snd) testData
        size = assert (n P.== n2) $ Z :. n
        vis = fromList size $ P.fst testData
        uvw = fromList size $ P.snd testData
        src      = undefined
        theta    = 2*0.05
        lam      = 18000
        kwargs   = noArgs {wstep= Just 2000, qpx= Just 1, npixFF= Just 256, npixKern= Just 31}
        otargs   = noOtherArgs
        (d,p, _) = CPU.run $ do_imaging theta lam uvw src vis w_cache_imaging kwargs otargs
        convKern = CPU.run $ w_kernel (constant theta) (constant 6000) kwargs
    

        fst i j = indexArray d (Z :. i :. j)
        maxi = P.maximum (toList d)
    --P.writeFile "data/mirror_uvw.csv" (makeVFile (showTriple $ printf "%e") (CPU.run muvw))
    --P.writeFile "data/result.csv" (makeMFile (printf "%e") d)
    --P.writeFile "data/result_p.csv" (makeMFile (printf "%e") p)
    P.putStrLn (P.show ((convKern)))

instance (P.Ord a) => P.Ord (Complex a)
    where
        (x1 :+ y1) <=  (x2 :+ y2) | x1 P.== x2 = y1 P.<= y2
                                  | P.otherwise = x1 P.<= x2

makeMFile :: (a -> String) -> Matrix a -> String
makeMFile showf mat = let
    (Z :. n :. m) = arrayShape  mat
    ls =  P.map showf $ toList mat
    lss = toRows n [] ls :: [[String]]
    rows = P.map (intercalate ",") lss :: [String]

    toRows _ res [] = P.reverse res
    toRows n res ys = let (x, xs) = P.splitAt n ys
                      in toRows n (x : res) xs

    in P.unlines rows

makeVFile :: (a -> String) -> Vector a -> String
makeVFile showf v = let
    (Z :. n) = arrayShape v
    ls =  P.map showf $ toList v
    in P.unlines ls

makeVFileC :: Vector (Complex Double) -> String
makeVFileC = makeVFile showC

makeMFileC :: Matrix (Complex Double) -> String
makeMFileC = makeMFile showC
    
showC :: Complex Double -> String
showC (x :+ y) | y P.>= 0 = '(' : printf "%e" x P.++ '+' : printf "%e" y P.++ "j)"
                | P.otherwise = '(' : printf "%e" x P.++ printf "%e" y P.++ "j)"

showTriple :: (a -> String) -> (a,a,a) -> String
showTriple print (x,y,z) = (intercalate ",") $ print x : print y : print z : []
-------------------------
-- Test input reading
testData0 :: IO ([Complex Double], [(Double, Double, Double)])
testData0 = do
    visData <- P.readFile "vis.csv"
    uvwData <- P.readFile "uvw.csv"
    let vislines = P.lines visData
        visdat   = P.map parseComplex vislines
        uvwlines = P.lines uvwData
        uvwdat   = P.map parseTuple3 uvwlines
    return (visdat, uvwdat)

parseComplex :: String -> Complex Double
parseComplex s  =
    case P.head s of
        ' ' -> parseComplex (P.tail s)
        '(' -> parseComplex (P.tail s)
        _   -> let (x, rest) = readDouble s
                   (y, _)    = readDouble rest
               in P.read x :+ P.read y
    where
        readDouble s@('-':_) = P.splitAt 25 s
        readDouble ('+':s)   = P.splitAt 24 s
        readDouble s         = P.splitAt 24 s

parseTuple3 :: P.Read a => String -> (a,a,a)
parseTuple3 s = 
    let valid c = c P./= ',' P.&& P.not (isSpace c)
        (x, _ : rest0) = P.span valid s
        (y, _ : rest1) = P.span valid rest0
        (z, _)         = P.span valid rest1
    in (P.read x, P.read y, P.read z)

testing :: Acc (Vector Int)
testing = use $ fromList (Z :. 3) [0..]

testingM :: Acc (Matrix Int)
testingM = use $ fromList (Z :. 3 :. 3) [0..]

reptesting :: Acc (Array DIM3 Int)
reptesting = replicate (constant (Z :. All :. (3::Int) :. (3:: Int))) $ testing

testing2 :: Exp Int -> Exp Int -> Exp Int -> (Exp Int)
testing2 a b c = reptesting ! index3 a b c

runExp :: Elt e => Exp e -> e
runExp e = indexArray (CPU.run (unit e)) Z

testing0 :: Vector (Int)
testing0 = fromList (Z :. 5) [0..]

testing3 :: Int
testing3 = indexArray testing0 (Z :. 4)

testing4 :: Matrix Int
testing4 = fromList (Z :. 5 :. 2) [0..]

testing5 :: Int -> Int -> Int
testing5 x y = runExp $ (use testing4) ! index2 (constant y) (constant x)

testing6 :: Elt e => Acc (Matrix e) -> Acc (Matrix e)
testing6 guv = let
    sh = shape guv
    Z :. n :. _ = unlift sh :: Z :. Exp Int :. Exp Int
    indexer (unlift -> Z :. y :. x) = if x== 0 || y == 0 
        then index2 y x 
        else index2 (n - y) (n - x) 

    oddReshape = (reverseOn _2 . reverseOn _1) guv
    evenReshape = backpermute sh indexer guv
    in if even n then evenReshape else oddReshape

ffttest0 :: Int -> (Vector Double, Vector Double, Vector Double, Vector Double)
ffttest0 n = let
    halfn = fromIntegral $ (n - 1) `div` 2
    minstart = case P.even n of
        True -> -halfn - 1
        False -> -halfn
    lst = [0..halfn] P.++ [minstart..(-1)]
    freqs = fromList (Z :. n) lst

    res  =  shift1D (use freqs)
    res2 = ishift1D (use freqs)
    res3 = ishift1D (res)
    --res3 = shiftTest lst
    in (freqs, CPU.run res,CPU.run res2, CPU.run res3)

ffttest1 :: Int -> Int -> (Matrix Double, Matrix Double, Matrix Double, Matrix Double)
ffttest1 n0 n1 = let
    n = n0 P.* n1
    halfn = fromIntegral $ (n - 1) `div` 2
    minstart = case P.even n of
        True -> -halfn - 1
        False -> -halfn
    lst = [0..halfn] P.++ [minstart..(-1)]
    freqs = fromList (Z :. n0 :. n1) lst

    res  =  shift2D (use freqs)
    res2 = ishift2D (use freqs)
    res3 = ishift2D (res)
    --res3 = shiftTest lst
    in (freqs, CPU.run res,CPU.run res2, CPU.run res3)

ffttest1_ :: Int -> Int -> Matrix Double
ffttest1_ n0 n1 = let
    n = n0 P.* n1
    halfn = fromIntegral $ (n - 1) `div` 2
    minstart = case P.even n of
        True -> -halfn - 1
        False -> -halfn
    lst = [0..halfn] P.++ [minstart..(-1)]
    in fromList (Z :. n0 :. n1) lst

ffttest2 :: Int -> Int -> Int -> (Array DIM3 Double, Array DIM3 Double, Array DIM3 Double, Array DIM3 Double)
ffttest2 n0 n1 n2 = let
    n = n0 P.* n1 P.* n2
    halfn = fromIntegral $ (n - 1) `div` 2
    minstart = case P.even n of
        True -> -halfn - 1
        False -> -halfn
    lst = [0..halfn] P.++ [minstart..(-1)]
    freqs = fromList (Z :. n0 :. n1 :. n2) lst

    res  =  shift3D (use freqs)
    res2 = ishift3D (use freqs)
    res3 = ishift3D (res)
    --res3 = shiftTest lst
    in (freqs, CPU.run res,CPU.run res2, CPU.run res3)

shiftTest :: [Double] -> [(Double, Int)]
shiftTest xs = let
    temp = addindex 0 xs
    addindex n [] = []
    addindex n (x:xs) = (x, n) : addindex (n+1) xs

    w = P.length xs
    mw      = w `P.div` 2
    mw2     = (w + 1) `P.div` 2
    indexer (v,x) = case x P.< mw2 of
        True  -> (v, x + mw)
        False -> (v, x -mw2)
    in P.map indexer temp


padTest :: Acc (Matrix Int)
padTest = padder testingM diml diml 0
        where
            dim = (1,1) ::  (Exp Int, Exp Int)
            diml = lift dim

testingCon :: Acc (Matrix Int)
testingCon = concatOn _2 m1 m2
    where
        m1 = use $ fromList (Z :. 1 :. 5) [0..]
        m2 = use $ fromList (Z :. 1 :. 5) [5..]


thekernelsTest :: Acc (Matrix Int)
thekernelsTest = let cpusteps = 16
                     makeWKernel i = use $ fromList (Z :. 1 :. 2) [i, i]
    in myfor (cpusteps - 1) (\i old -> concatOn _2 (makeWKernel i) old) (makeWKernel (cpusteps - 1))