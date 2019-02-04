{-# language FlexibleContexts    #-}
{-# language RebindableSyntax    #-}
{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}

module Main where

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A
import Data.Array.Accelerate.Math.DFT.Centre              as A
import Data.Array.Accelerate.Math.FFT                     as A

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I
import Data.Array.Accelerate.Debug                  as A

import Gridding

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..))
import Data.Char (isSpace)
import Text.Printf
import Data.List (intercalate)

import Debug.Trace
import System.IO.Unsafe (unsafePerformIO)
import Control.Exception (assert)

import Control.Lens as L (_1, _2, (^.) )

main :: IO ()
main = testConv

setflags :: IO ()
setflags = do
    --setFlag dump_cc
    --setFlag dump_sched
    setFlag dump_phases

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
        thecode = do_imaging theta lam uvw src vis simple_imaging kwargs otargs
        (d,p, _) = CPU.run $ thecode
        (muvw, mvis) = unzip $ mirror_uvw (use uvw) (use vis)
        (dmax, pmax) = CPU.run $ dpmax thecode
    --P.writeFile "data/mirror_uvw.csv" (makeVFile (showTriple $ printf "%e") (CPU.run muvw))
    --P.writeFile "data/result.csv" (makeMFile (printf "%e") d)
    --P.writeFile "data/result_p.csv" (makeMFile (printf "%e") p)
    --P.writeFile ("data/generated_simple.hs") (P.show thecode)
    P.putStrLn (P.show dmax)

testWCache :: IO ()
testWCache = do
    testData <- testData0
    let n   = (P.length . P.fst) testData
        n2   = (P.length . P.snd) testData
        size = assert (n P.== n2) $ Z :. n
        vis = fromList size $ P.fst testData
        uvw = fromList size $ P.snd testData
        src      = undefined
        theta    = 2*0.005
        lam      = 18000
        qpx      = 2
        w        = 10000
        wstep    = constant w
        npixFF   = 256
        npixKern = 31
        kwargs   = noArgs {wstep= Just w, qpx= Just qpx, npixFF= Just npixFF, npixKern= Just npixKern}
        otargs   = noOtherArgs
        --convKern = CPU.run $ w_kernel (constant theta) (constant (fromIntegral w)) kwargs
        --convKern' = arrayReshape (Z :. (arraySize . arrayShape ) convKern ) convKern
        dimg = do_imaging theta lam uvw src vis w_cache_imaging kwargs otargs
        (d,p, pmax2) = CPU.run $ dimg
        (dmax, pmax) = CPU.run $ dpmax dimg

        (l, m) = unlift $ kernel_coordinates (constant 256) (constant theta) kwargs :: (Acc (Matrix Double), Acc (Matrix Double))
        kern = w_kernel_function l m (constant (fromIntegral w))
        padff = pad_mid kern (constant npixFF * constant qpx)
        af = ifft padff

        (_,_,ww) = unzip3 (use uvw)
        roundedw = map (\w' -> wstep * round ( w' / A.fromIntegral wstep)) ww
        maxww = the $ maximum ww 
        maxw = the $ maximum roundedw :: Exp Int
        minw = the $ minimum roundedw :: Exp Int
        steps = ((maxw - minw) `div` wstep ) + 1 :: Exp Int
        (cpumin, cpumax, cpusteps, cpumaxww) = CPU.run (unit $ lift (minw, maxw, steps, maxww)) `indexArray` Z
        wbins = map (\rw' -> (rw' - constant cpumin) `div` wstep) roundedw
        

    --P.writeFile "data/mirror_uvw.csv" (makeVFile (showTriple $ printf "%e") (CPU.run muvw))

    P.writeFile "data/result.csv" (makeMFile (printf "%e") d)
    P.writeFile "data/result_p.csv" (makeMFile (printf "%e") p)
    P.writeFile "data/roundedw.csv" (makeVFile (P.show) (CPU.run roundedw))
    --P.putStrLn (P.show ((dmax, pmax)))
    --P.writeFile "data/result_conv.csv" (makeVFile (showC) convKern')
    --P.writeFile "data/result_kern.csv" (makeMFile (showC) kern)
    --P.writeFile "data/result_af.csv" (makeMFile (showC) (CPU.run af))
    --P.putStrLn (P.show ((convKern `indexArray` (Z :. 0 :. 0 :. 0 :. 0) )))
    P.putStrLn (P.show (cpumin, cpumax, cpusteps, cpumaxww) )
    P.putStrLn (P.show pmax2)
    --P.putStrLn (P.show (CPU.run wbins))

testConv :: IO ()
testConv = do
    testData <- testData0
    let n   = (P.length . P.fst) testData
        n2   = (P.length . P.snd) testData
        size = assert (n P.== n2) $ Z :. n
        vis = fromList size $ P.fst testData
        uvw = fromList size $ P.snd testData
        src      = undefined
        theta    = 2*0.05
        lam      = 18000
        lamf     = fromIntegral lam
        qpx      = 2
        w        = 1000
        npixFF   = 256
        npixKern = 31
        ne = constant $ P.round(theta * lamf) :: Exp Int
        convKern = CPU.run $ w_kernel (constant theta) (constant (fromIntegral w)) kwargs
        convKern' = CPU.run $ flatten (use convKern)
        kwargs   = noArgs {wstep= Just w, qpx= Just qpx, npixFF= Just npixFF, npixKern= Just npixKern}
        otargs   = noOtherArgs {convolutionKernel = Just convKern}
        dimg = do_imaging theta lam uvw src vis conv_imaging kwargs otargs 
        (d,p, _) = CPU.run dimg
        (dmax, pmax) = CPU.run $ dpmax dimg
        
    --P.writeFile "data/mirror_uvw.csv" (makeVFile (showTriple $ printf "%e") (CPU.run muvw))
    --P.writeFile "data/result.csv" (makeMFile (printf "%e") d)
    --P.writeFile "data/result_coords.csv" (makeVFile (showQuad P.show) (CPU.run coords))
    --P.writeFile "data/resultConv_kern.csv" (makeVFile (showC) convKern')
    --P.writeFile "data/result_p.csv" (makeMFile (printf "%e") p)
    --P.writeFile "data/resultConv.csv" (makeMFile (printf "%e") d)
    P.putStrLn (P.show ((dmax, pmax)))

instance (P.Ord a) => P.Ord (Complex a)
    where
        (x1 :+ y1) <=  (x2 :+ y2) | x1 P.== x2 = y1 P.<= y2
                                  | P.otherwise = x1 P.<= x2

makeMFile :: Elt a => (a -> String) -> Matrix a -> String
makeMFile showf mat = let
    (Z :. n :. m) = arrayShape  mat
    ls =  P.map showf $ toList mat
    lss = toRows n [] ls :: [[String]]
    rows = P.map (intercalate ",") lss :: [String]

    toRows _ res [] = P.reverse res
    toRows n res ys = let (x, xs) = P.splitAt n ys
                      in toRows n (x : res) xs

    in P.unlines rows

makeVFile :: Elt a => (a -> String) -> Vector a -> String
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

showQuad :: (a -> String) -> (a,a,a,a) -> String
showQuad print (w,x,y,z) = (intercalate ",") $ print w : print x : print y : print z : []
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

testFixbounds :: Acc (Matrix (Complex Double))-- Acc (Vector (Int, Int, Int))
testFixbounds =
    let
        v1 ::  Acc (Vector (Int, Int, Complex Double))
        v1 = use $ fromList (Z :. 10) [((2 * x)  `P.mod` 5 ,(3 * x + 1) `P.mod` 5, ( (fromIntegral $ x+5) :+ 1.0)) | x <- [0..] ]
        width = constant 5 :: Exp Int
        height = constant 5 :: Exp Int
        cc0 = (lift (constant 0.0 :+ constant 0.0))
        a = fill (index2 height width) cc0

        fullrepli i j origin =
            let
                offsety = constant $ i :: Exp Int
                offsetx = constant $ j :: Exp Int
                f = fixoutofboundsOLD offsetx offsety width height cc0
                checkbounds = map f v1
                (x, y, source) = unzip3 $ checkbounds :: (Acc (Vector Int), Acc (Vector Int), Acc (Vector (Complex Double)))
                indexer id =
                    let y' = y ! id + offsety
                        x' = x ! id + offsetx
                    in lift (Z :. y' :. x')
            in permute (+) origin indexer source
    in fullrepli 0 1 $fullrepli 1 0 $ fullrepli 1 1 $ fullrepli 0 0 a

testFixbounds2 :: Acc (Matrix (Complex Double))-- Acc (Vector (Int, Int, Int))
testFixbounds2 =
    let
        v1 ::  Acc (Vector (Int, Int, Int, Int, Complex Double))
        v1 = use $ fromList (Z :. 10)
            [((2 * x)  `P.mod` 5
            , x `P.mod` 2
            ,(3 * x + 2) `P.mod` 5
            , x `P.mod` 2
            , ( (fromIntegral $ x+5) :+ 1.0))
                | x <- [0..] ]
        width = constant 5 :: Exp Int
        height = constant 5 :: Exp Int
        cc0 = (lift (constant 0.0 :+ constant 0.0))
        a = fill (index2 height width) cc0

        fixer = fixoutofbounds width height cc0

        gcf = use $ fromList (Z :. 2 :. 2 :. 2 :. 2) [ x :+ x | x <- [0..] ] :: Acc (Array DIM4 (Complex Double) )

        getComplex :: Exp DIM3 -> Exp (Int, Int, Int, Int, Complex Double) -> Exp (Int, Int, Complex Double)
        getComplex
            (unlift . unindex3 -> (_, i, j) ::(Exp Int,Exp Int,Exp Int))
            (unlift -> (x, xf, y, yf, vis)::(Exp Int,Exp Int,Exp Int,Exp Int,Exp (Complex Double))) =
            lift ( x + j
                 , y + i
                 , vis * gcf ! lift (Z :. yf :. xf :. i :. j) )


        totv1 = replicate (constant (Z :. All :. (2 :: Int) :. (2 :: Int))) v1
        withkern = imap getComplex totv1
        (x,y, val) = unzip3 $ map fixer withkern

        indexer id =
            let y' = y ! id
                x' = x ! id
            in index2 y' x'
    in permute (+) a indexer val


fixoutofboundsOLD ::Elt a => Exp Int -> Exp Int -> Exp Int -> Exp Int -> Exp a -> Exp (Int, Int, a) -> Exp (Int, Int, a)
fixoutofboundsOLD offsetx offsety maxx maxy def
                old@(unlift -> (x', y', _) ::(Exp Int,Exp Int,Exp a)) =
    let idx = x' + offsetx
        idy = y' + offsety
        outofbounds = idx < c0 || idy < c0 || idx >= maxx || idy >= maxy
        -- If out of bounds, set all values to zero (after the offset alters it)
        newx = - offsetx :: Exp Int
        newy = - offsety :: Exp Int
        newv = def
        new = lift (newx, newy, newv)
    in if outofbounds then new else old

dpmax :: Acc (Matrix Double,Matrix Double, Scalar Double) ->  Acc (Scalar Double, Scalar Double)
dpmax (unlift -> (d, _, pmax) :: (Acc (Matrix Double),Acc (Matrix Double),Acc (Scalar Double)) ) =
    let dmax = maximum . flatten $ d
    in lift (dmax, pmax) 

testConcat :: Acc (Array DIM2 Int)
testConcat =
    let mylist i = let myi = use $ fromList Z [i]
                   in replicate (constant (Any :. (1::Int) :. (4::Int))) myi
        cpusteps = 3
    in myfor (cpusteps - 1) (\i old -> concatOn _2 (mylist i) old) (mylist (cpusteps - 1))