{-# LANGUAGE CPP #-}
module Main where

import ImageDataset
import FFT
import System.Directory
import Control.Exception
import System.IO.Error hiding (catch)
import System.Environment
import qualified Data.Array.Accelerate.LLVM.Native      as CPU
import qualified Data.Array.Accelerate.Interpreter      as I
import Data.Array.Accelerate.Debug

#ifdef ACCELERATE_LLVM_PTX_BACKEND
import qualified Data.Array.Accelerate.LLVM.PTX         as GPU
#endif

data Args = Args { n   :: Maybe Int
                 , runner :: Run
                 , input :: String
                 , out :: Maybe String
                 , flags :: [String]
                 , chunks :: Maybe Int
                 , fftf   :: FFT     }

defArgs = Args (Just 1) CPU "data" Nothing [] Nothing Adhoc

main :: IO ()
main = do
    --remover "result.h5"
    args <- getArgs
    let parsedArgs = parser defArgs args
        n_         = n parsedArgs
        backend    = runner parsedArgs
        dir        = input parsedArgs
        wkern      = dir ++ "/" ++ "SKA1_Low_wkern2.h5"
        akern      = dir ++ "/" ++ "SKA1_Low_akern3.h5"
        visdata    = dir ++ "/" ++ "SKA1_Low_quick.h5"
        outFile    = out parsedArgs
        fs         = flags parsedArgs
        chunks_    = chunks parsedArgs
        fft_       = fftf parsedArgs
        runner1 :: Runners
        runner1    = case backend of
            CPU   -> CPU.run1
            Inter -> I.run1
#ifdef ACCELERATE_LLVM_PTX_BACKEND
            GPU   -> GPU.run1
#else
            GPU -> error "GPU implementation not supported, build with flag \"llvm-gpu\""
#endif
    case chunks_ of
        Nothing -> return ()
        Just n  -> do putStrLn ("Chunks of " ++ show n); setEnv "ACCELERATE_FLAGS" ("-chunk-size=" ++ show n)
    
    myfft<- case fft_ of
              Adhoc   -> do putStrLn "Adhoc fft"; return fft2DAdhoc
              Old     -> do putStrLn "Old fft"; return fft2DOld
              Foreign -> do putStrLn "Foreign fft";return fft2DFor
        
    case backend of
        CPU   -> putStrLn "CPU backend"
        Inter -> putStrLn "Interpreter backend"
        GPU   -> putStrLn "GPU backend"
    mapM_ processFlags fs
    
    case outFile of
        Nothing -> return ()
        Just fn  -> remover fn
    fourier <- aw_gridding runner1 myfft backend fft_ wkern akern visdata n_ outFile
    putStrLn (show fourier)


remover :: FilePath -> IO ()
remover fileName = removeFile fileName `catch` handle
    where handle e
            | isDoesNotExistError e = return ()
            | otherwise = throwIO e

parser :: Args -> [String] -> Args
parser args []            = args
parser args ("-debug":xs) = parser args{runner = Inter} xs
parser args ("-g":xs)     = parser args{runner = GPU} xs
parser args ("-G":xs)     = parser args{runner = GPU} xs
parser args ("-gpu":xs)   = parser args{runner = GPU} xs
parser args ("-GPU":xs)   = parser args{runner = GPU} xs
parser args ("-Gpu":xs)   = parser args{runner = GPU} xs
parser args ("-old":xs)   = parser args{fftf = Old} xs
parser args ("-for":xs)   = parser args{fftf = Foreign} xs
parser args ("-foreign":xs)
                          = parser args{fftf = Foreign} xs
parser args ("-adhoc":xs) = parser args{fftf = Adhoc} xs
parser args ("-n":n:xs)   = parser args{n = Just $ read n} xs
parser args ("-all":xs)   = parser args{n = Nothing} xs
parser args ("-i":inp:xs) = parser args{input = inp} xs
parser args ("-o":out:xs) = parser args{out = Just out} xs
parser args ("-chunks":n:xs) = parser args{chunks = Just $ read n} xs
parser args@Args{flags = fs} (('-':'d':f):xs) = parser args{flags=f:fs} xs
parser args xs            = error $ "Error while parsing" ++ show xs

processFlags :: String -> IO()
processFlags "verbose" = setFlag verbose
processFlags "dump_phases" = setFlag dump_phases
processFlags "dump_sharing" = setFlag dump_sharing
--processFlags "dump_fusion" = setFlag dump_fusion
processFlags "dump_simpl_stats" = setFlag dump_simpl_stats
processFlags "dump_simpl_iterations" = setFlag dump_simpl_iterations
processFlags "dump_vectorisation" = setFlag dump_vectorisation
processFlags "dump_dot" = setFlag dump_dot
processFlags "dump_simpl_dot" = setFlag dump_simpl_dot
processFlags "dump_gc" = setFlag dump_gc
processFlags "dump_gc_stats" = setFlag dump_gc_stats
processFlags "dump_cc" = setFlag dump_cc
processFlags "dump_ld" = setFlag dump_ld
processFlags "dump_asm" = setFlag dump_asm
processFlags "dump_exec" = setFlag dump_exec
processFlags "dump_sched" = setFlag dump_sched
processFlags xs = putStrLn $ "Warning: Flag \'" ++ xs ++ "\' isn't recognized"
