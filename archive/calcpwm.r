tmp <- matrix(c(1L,2L,2L,0L, 2L,1L,2L,0L, 3L,0L,1L,1L, 0L,5L,0L,0L, 5L,0L,0L,0L, 0L,0L,4L,1L, 0L,1L,4L,0L, 0L,0L,0L,5L, 0L,0L,5L,0L, 0L,1L,2L,2L, 0L,2L,0L,3L, 1L,0L,3L,1L), nrow=4, byrow=FALSE, dimnames=list(c("A", "C", "G", "T")))
tseq <- "ttctggttatggcttctttcacctcaacctctccaagacccttattttcctcacctctcatccaaatccacattctacttgtctctcattagaaaactaacaggtatccaaagaacaacattaaTTTCTATACTTTTGAAGCTAGAATATCCTGAATATAAACTTTGAATTGCCATATTAATATTTATTTTCAACAGACTATGgggaaaccagaatgctaagtcaacggaatagaacagattttgacaagatcatagaggaaaacttccccaaactcaggaaagacatacctatacccctcaacaagaaggaaaccatgcctggtaatggaaacttagctagctacccagaactaatgatatcatggatcttggggcagaaatttttctaaggaagagcacaccagttgattagccaatatgaaaatgtcagccctgaaaatgtacacatataagacagactttgcaggttatatttatgtatttagATACACtgttgggagccattaaggcaacgctattgtcctgatctctgaattgggcctctccccccaagaagaaaagagggtcaaaagcgggccaccgacacaccgctctgagaacaaccacagattgttccagccctaagtcagcgccggatgtcctgacctcaagataccctgatactgccaagttcctgcttccccatgtagctgccaaaagaaagaatttctattgccccccctcccataagtacttctccttttgcttgtatttccccatctctccactgataagtatctgtacttccccgtttccttgtgcaattctactgccccatccttcctgctgaaaggcacttccccttttgcttgtgtatttaaaccttgggcctggctaatacatttgggggtcttgatacaacttcagaatgctttctgtgtcgttattcgtacaagaccctcgtctctctctaccccccatttggttattaggaggaggtcccccgagatcctcgaataactggacctgctggacgggtcaagtggcgcccgacgtggggcacgaggtacgactaccctccggacaacggattaaggattaaaaggtaccgcgcagacagagaaacagtgggacagtccaccgacacttcactcgagtgtcgagtacggacgtcgggaggtaaacctaggcattcagttgttgtcccattaaccatggggcatttagtgaagatataaagggttctcttaaggagagagggataagagttttaaaaaaaggggaatttaataatatatttttgttttgtggacaagagatgtccatggctggtattaaatggttctgaaacttagaataaagtaggcaaagaaatcaatagtttaataaaacaaagagaagatgtccatgagtccttttttagttattggggaatcatcagagatctccttgagcaggcagaaaagggcgcaaaaggtgctcgcctcctagcccttactgaagattttttagaaagctcccaatcgccatctcacaaaagcaaaactaagctgccagcgcacagtatctctattagaaggtttagcagaaggctacttggtacccaataaatgaaataaagttattcagtctgtcctcactagaggccaatatttaacttaaaagtcagagtttGtggaagacagcctcctggactgctgataaaatctgtggttaaggcaatttcgctgcagacagaaaacaattagggctctcccctattgtcctcgcccagacagcacaagcggccctgagagcttggcgcgctgtctctgcagcataggggcactaactatgcccctaaccaagattatacaagggccccaagaatcctatgcacagttcgttgacagattgcaggaagcagctgagagaattttgggacctgaggaaagtgagggtctgttagttcgacaacttgctcttgaaaacgccaattcagcatgtagggctgccctgagaggtaaaaccaagagttta"
res <- matchPWM(PWM(tmp), tseq, min.score="0%")
res <- PWMscoreStartingAt(PWM(tmp), tseq, 1:100)
print(res)

inputFile <- "/home/path/to/myFile.txt"
con  <- file(inputFile, open = "r")

dataList <- list()
ecdfList <- list()

while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    myVector <- (strsplit(oneLine, " "))
    myVector <- list(as.numeric(myVector[[1]]))
    dataList <- c(dataList,myVector)

    myEcdf <- ecdf(myVector[[1]])
    ecdfList <- c(ecdfList,myEcdf)

  } 

close(con)