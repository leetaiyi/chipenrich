context('Test file reading logic')

# Wrap the reading logic from chipenrich(...) into a function to ease testing
reading_logic = function(peaks) {
	if (class(peaks) == "data.frame") {
		peakobj = load_peaks(peaks);
	} else if (class(peaks) == "character") {
    if (stringr::str_sub(peaks,-4,-1) == ".gff" || stringr::str_sub(peaks,-5,-1) == '.gff3' || stringr::str_sub(peaks,-7,-1) == ".gff.gz" || stringr::str_sub(peaks,-8,-1) == '.gff3.gz') {
      message("Reading peaks file: ",peaks);
      peakobj = read_bedgff(peaks);
    } else {
      message("Reading peaks file: ",peaks);
      peakobj = read_bed(peaks);
    }
	}
  return(peakobj)
}

peak_path = system.file('extdata', 'test.gff3', package='chipenrich')
peaks = reading_logic(peak_path)
test_that('.gff3 gives correct number of records after reduce peaks',{
    expect_equal(sum(sapply(peaks,length)), 19)
})

peak_path = system.file('extdata', 'test.gff3.gz', package='chipenrich')
peaks = reading_logic(peak_path)
test_that('.gff3.gz gives correct number of records after reduce peaks',{
    expect_equal(sum(sapply(peaks,length)), 19)
})

peak_path = system.file('extdata', 'test_header.bed', package='chipenrich')
peaks = reading_logic(peak_path)
test_that('.bed with (2 overlapping peaks) header gives correct number of records after reduce peaks',{
    expect_equal(sum(sapply(peaks,length)), 9)
})

peak_path = system.file('extdata', 'test.broadPeak', package='chipenrich')
peaks = reading_logic(peak_path)
test_that('.broadPeak (1 duplicate line and 2 overlapping peaks) gives correct number of records after reduce peaks',{
    expect_equal(sum(sapply(peaks,length)), 25)
})

peak_path = system.file('extdata', 'test.narrowPeak', package='chipenrich')
peaks = reading_logic(peak_path)
test_that('.narrowPeak gives correct number of records after reduce peaks',{
    expect_equal(sum(sapply(peaks,length)), 14)
})

peak_path = system.file('extdata', 'test.broadPeak.gz', package='chipenrich')
peaks = reading_logic(peak_path)
test_that('.broadPeak.gz gives correct number of records after reduce peaks',{
    expect_equal(sum(sapply(peaks,length)), 17)
})

peak_path = system.file('extdata', 'test.narrowPeak.gz', package='chipenrich')
peaks = reading_logic(peak_path)
test_that('.narrowPeak.gz gives correct number of records after reduce peaks',{
    expect_equal(sum(sapply(peaks,length)), 11)
})
