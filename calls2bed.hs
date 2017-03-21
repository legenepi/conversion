import System.Environment
import System.IO
import qualified Data.ByteString.Char8 as C
import Data.ByteString.Lazy.Builder
import Data.Bits
import GHC.Word
import Data.Monoid

homozygote1 :: Word8
homozygote1 = 0

heterozygote :: Word8
heterozygote = 2

homozygote2 :: Word8
homozygote2 = 3

missing :: Word8
missing = 1

magicNo :: Word16
magicNo = 27675

mode :: Word8
mode = 1

procHead :: Handle -> Handle -> IO ()
procHead hC hF = do
    line <- C.hGetLine hC
    if C.head line == '#'
        then procHead hC hF
        else C.hPut hF $ getIds line where
            getIds = C.unlines . tail . C.words

writeBedHeader :: Handle -> IO ()
writeBedHeader hB = hPutBuilder hB header where
    header = word16BE magicNo <> word8 mode

procLines :: Handle -> Handle -> Handle -> IO ()
procLines hC hB hI = do
    eof <- hIsEOF hC
    if eof
        then return ()
        else do
            line <- C.hGetLine hC
            procLine hB hI $ C.words line 
            procLines hC hB hI

procLine :: Handle -> Handle -> [C.ByteString] -> IO ()
procLine hB hI (x:xs) = do
    C.hPutStrLn hI x
    hPutBuilder hB $ procGenos mempty xs

procGenos :: Builder -> [C.ByteString] -> Builder
procGenos b [] = b
procGenos b g = procGenos bytes rest where
        (geno, rest) = splitAt 4 g
        bytes = b <> (bits2bytes . geno2bits $ geno)

        geno2bits :: [C.ByteString] -> [Word8]
        geno2bits x = zipWith shiftL y [0, 2, 4, 6] where
            y = map int2bits x
        
        bits2bytes :: [Word8] -> Builder
        bits2bytes bits = word8 $ foldl1 (+) bits

        int2bits :: C.ByteString -> Word8
        int2bits is
            | C.elem '0' is = homozygote1
            | C.elem '1' is = heterozygote
            | C.elem '2' is = homozygote2
            | otherwise = missing

main :: IO ()
main = do
    args <- getArgs
    hCalls <- openFile (args !! 0) ReadMode
    hBed <- openFile "out.bed" WriteMode
    hBim <- openFile "out.bim" WriteMode
    hFam <- openFile "out.fam" WriteMode
    procHead hCalls hFam
    hClose hFam
    writeBedHeader hBed
    procLines hCalls hBed hBim
    hClose hCalls
    hClose hBed
    hClose hBim
