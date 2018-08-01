package main

import (
	//	"bufio"
	"flag"
	"fmt"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/seqyuan/go-intervals/intervalset"
	"io"
	"log"
	"os"
	"strings"
	"math"
	"strconv"
    "path/filepath"
)

func usage() {
	//Usage function
    fmt.Println("Usage : mapped_2hic_fragments")
    fmt.Println("-f/--fragmentFile <Restriction fragment file GFF3>")
    fmt.Println("-r/--mappedReadsFile <BAM/SAM file of mapped reads>")
    fmt.Println("[-o/--outputDir] <Output directory. Default is current directory>")
    fmt.Println("[-s/--shortestInsertSize] <Shortest insert size of mapped reads to consider>")
    fmt.Println("[-l/--longestInsertSize] <Longest insert size of mapped reads to consider>")
    fmt.Println("[-t/--shortestFragmentLength] <Shortest restriction fragment length to consider>")
    fmt.Println("[-m/--longestFragmentLength] <Longest restriction fragment length to consider>")
    fmt.Println("[-d/--minCisDist] <Minimum distance between intrachromosomal contact to consider>")
    fmt.Println("[-g/--gtag] <Genotype tag. If specified, this tag will be reported in the valid pairs output for allele specific classification>")
    fmt.Println("[-a/--all] <Write all additional output files, with information about the discarded reads (self-circle, dangling end, etc.)>")
    fmt.Println("[-S/--sam] <Output an additional SAM file with flag 'CT' for pairs classification>")
    fmt.Println("[-v/--verbose] <Verbose>")
    fmt.Println("[-h/--help] <Help>")
}


func get_read_strand(r *sam.Record) (strand string) {
	/*
	   Conversion of read position to naive strand representation

	   Parameters
	   ----------
	   read : list
	       list of aligned reads
	*/
	strand = "+"
	if r.Strand() == -1 {
		strand = "-"
	}
	return
}

func isIntraChrom(r1 *sam.Record, r2 *sam.Record) bool {
	/*
	   Return true is the reads pair is intrachromosomal

	   read1 : [AlignedRead]
	   read2 : [AlignedRead]

	*/

	if r1.Ref.Name() == r2.Ref.Name() {
		return true
	}
	return false
}

def get_read_pos(read *sam.Record)(Pos int){
    /*
    Return the read position (zero-based) used for the intersection with
    the restriction fragment

    The 5' end is not a good choice for the reverse reads (which contain part
    of the restriction site, and thus overlap the next restriction fragment)
    Using the left-most position (5' for forward, 3' for reverse) or the
    middle of the read should work but the middle of the reads might be more
    safe

    Parameters
    -----------
    read : list
        list of aligned reads
    */
    pos = read.Pos + int(read.Len()/2)
    return Pos
}

func get_cis_dist(r1 *sam.Record, r2 *sam.Record)dist int{
     /*
     Calculte the contact distance between two intrachromosomal reads
     """
     # Get oriented reads
     ##r1, r2 = get_ordered_reads(read1, read2)
     */
     dist = -1
     if strings.Contains(r1.Flags.String(), "u") || strings.Contains(r2.Flags.String(), "u"){         
         // Contact distances can be calculated for intrachromosomal reads only
     	return
     }

     if isIntraChrom(r1, r2){
        r1pos := get_read_pos(r1)
        r2pos := get_read_pos(r2)
        dist = math.Abs(r1pos - r2pos)
     }
     return
}

func get_read_start(read *sam.Record)pos int{ 
    //Return the 5' end of the read
    if read.Strand() == -1{
        pos = read.Pos + read.Len()
    }
    else:
        pos = read.Pos
    return
}

func get_ordered_reads(read1 *sam.Record, read2 *sam.Record)(r1 *sam.Record, r2 *sam.Record){
    /*
    Reorient reads

    The sequencing is usually not oriented. Reorient the reads so that r1 is
    always before r2.
    Sequencing is always performed from 5' to 3' end
    So in unstranded case, we can have

    1              2
    --->           --->
    ==========  or =========
         <----          <---
             2             1

    Reordering the reads allow to always be in the first case
    */
    if read1.Ref.Name() == read2.Ref.Name() {
        if get_read_pos(read1) < get_read_pos(read2):
            r1 = read1
            r2 = read2
        else:
            r1 = read2
            r2 = read1
    }else{
        if read1.RefID() < read2.RefID(){
            r1 = read1
            r2 = read2
        }else{
            r1 = read2
            r2 = read1
        }
    }
                
    return r1, r2
}

func load_restriction_fragment(in_file string, minfragsize int, maxfragsize int, verbose bool)map[string]*Set {
    /*
       Read a BED file and store the intervals in a tree
       Intervals are zero-based objects. The output object is a hash table with
       one search tree per chromosome
    */
    resFrag := make(map[string]*Set)

    if verbose {
        fmt.Println("## Loading Restriction File Intervals '", in_file, "'...")
    }

    rw, err := os.Open(in_file)
    if err != nil {
        panic(err)
    }
    defer rw.Close()
    rb := bufio.NewReader(rw)
    nline := 0
    for {
        nline += 1
        line_byte, _, err := rb.ReadLine()
        if err == io.EOF {
            break
        }
        bedtab := bytes.Split(line_byte, []byte{'\t'})
        var start int
        var end int

        //BED files are zero-based as Intervals objects
        name := string(bedtab[3])
        chromosome := string(bedtab[0])
        start, err = strconv.Atoi(string(bedtab[1]))
        end, err = strconv.Atoi(string(bedtab[2]))
        start += 1
        end += 1
        fragl := end - start

        // Discard fragments outside the size range
        if minfragsize != 0 && fragl < minfragsize {
            fmt.Println("Warning : fragment ", name, " [", fragl, "] outside of range. Discarded")
            continue
        }
        if maxfragsize != 0 && fragl > maxfragsize {
            fmt.Println("Warning : fragment ", name, " [", fragl, "] outside of range. Discarded")
            continue
        }

        frag_span := &Span{
            name,
            start,
            end,
        }


        if _, ok := resFrag[chromosome]; ok {
            tree := resFrag[chromosome]
            tree.DangerInsert(frag_span)
        } else {
            tree := Empty()
            tree.DangerInsert(frag_span)

            //tree := NewSet(frag_span)
            resFrag[chromosome] = tree
            fmt.Println(tree)
        }
    }
    return resFrag
}

func get_overlapping_restriction_fragment(resFrag map[string]*Set, chrom string, read *sam.Record) Interval{
    /*
    Intersect a given read with the set of restriction fragments

    resFrag = the restriction fragments [hash]
    chrom = the chromosome to look at [character]
    read = the read to intersect [AlignedRead]
	*/
    // Get read position (middle or 5' end)
    pos := get_read_pos(read)
    read_pos := &Span{
		"readsname",
		pos,
		pos+1,
	}

    read_pos_span := NewSet([]Interval{read_pos})
    
    if _, ok := resFrag[chrom]; ok{
        // Overlap with the position of the read (zero-based)
        read_pos_span.Intersect(resFrag[chrom])
        overlap_frags_len := read_pos_span.AllIntervals()
        if overlap_frags_len > 1{
            fmt.Println("Warning : ", overlap_frags_len, " restriction fragments found for ", read.Ref.Name(), "- skipped")
            return nil
        }else if overlap_frags_len == 0{
            fmt.Println("Warning - no restriction fragments for ", read.Ref.Name() ," at ", chrom, ":", pos)
            return nil
        }else{
            return read_pos_span.AllIntervals()[0]
        }
    }else{
        fmt.Println("Warning - no restriction fragments for ", read.Ref.Name()," at ", chrom, ":", pos)
        return nil
    }
}

func are_contiguous_fragments(frag1, frag2 Interval, chr1, chr2 string)ret bool{
    //Compare fragment positions to check if they are contiguous
    ret = false
    var d int
    if chr1 == chr2{
        if frag1.MIN() < frag2.MIN()){
            d = frag2.MIN() - frag1.MAX()
        }else{
            d = frag1.MIN() - frag2.MAX()
        }
            
        if d == 0{
            ret = true
        }
    }
    return
}

func is_religation(read1, read2 *sam.Record, frag1, frag2 Interval)ret bool{
    //Reads are expected to map adjacent fragments
    //Check the orientation of reads -><-
    ret = false
    if are_contiguous_fragments(frag1, frag2, read1.Name, read2.Name){
        //r1, r2 = get_ordered_reads(read1, read2)
        //if get_read_strand(r1) == "+" and get_read_strand(r2) == "-":
        ret = true
    }
    return
}

func is_self_circle(read1, read2 *sam.Record)bool{
    /*
    Both reads are expected to be on the same restriction fragments
    Check the orientation of reads <-->

    read1 : [AlignedRead]
    read2 : [AlignedRead]
    */
    ret := false
    //Get oriented reads
    r1, r2 = get_ordered_reads(read1, read2)
    //1<- ->2 or 2<- ->1
    if get_read_strand(r1) == "-" and get_read_strand(r2) == "+"{
        ret = true
    }
    return ret
}

func is_dangling_end(read1, read2 *sam.Record)bool{
    /*
    Both reads are expected to be on the same restriction fragments
    Check the orientation of reads -><-

    read1 : [AlignedRead]
    read2 : [AlignedRead]
    */
    ret := false
    // Get oriented reads
    r1, r2 = get_ordered_reads(read1, read2)
    // 1-> <-2 or 2-> <-1
    if get_read_strand(r1) == "+" and get_read_strand(r2) == "-":
        ret = true
    return ret
}

func get_valid_orientation(read1, read2 *sam.Record)string{
    /*
    Both reads are expected to be on the different restriction fragments
    Check the orientation of reads ->-> / <-<- / -><- / <-->

    read1 : [AlignedRead]
    read2 : [AlignedRead]
	*/
    // Get oriented reads
    r1, r2 = get_ordered_reads(read1, read2)

    var direction string
    if get_read_strand(r1) == "+" and get_read_strand(r2) == "+"{
        direction = "FF"
    }else if get_read_strand(r1) == "-" and get_read_strand(r2) == "-"{
        direction = "RR"
    }else if get_read_strand(r1) == "+" and get_read_strand(r2) == "-"{
        direction = "FR"
    }else if get_read_strand(r1) == "-" and get_read_strand(r2) == "+"{
        direction = "RF"
    }

    return direction
}

def get_PE_fragment_size(read1, read2 *sam.Record, resFrag1, resFrag2 Interval, interactionType string)int{
    /*
    Calculte the size of the DNA fragment library
    read1 : [AlignedRead]
    read2 : [AlignedRead]
    resfrag1 = restrictin fragment overlapping the R1 read [interval]
    resfrag1 = restrictin fragment overlapping the R1 read [interval]
    interactionType : Type of interaction from get_interaction_type() [str]
	*/

    var fragmentsize int
    var rfrag1 Interval
    var rfrag2 Interval

    // Get oriented reads
    r1, r2 := get_ordered_reads(read1, read2)
    
    if !strings.Contains(r1.Flags.String(), "u") && !strings.Contains(r2.Flags.String(), "u"):
        if r1 == read2:
            rfrag1 = resFrag2
            rfrag2 = resFrag1
        else:
            rfrag1 = resFrag1
            rfrag2 = resFrag2

        // In this case use the read start !
        r1pos := get_read_start(r1)
        r2pos := get_read_start(r2)

        if interactionType == "DE" || interactionType == "RE"{
            fragmentsize = r2pos - r1pos
        }else if interactionType == "SC"{
            fragmentsize = (r1pos - rfrag1.MIN()) + (rfrag2.MAX() - r2pos)
        }else if interactionType == "VI"{
        	var dr1,dr2 int
            if get_read_strand(r1) == "+"{
                dr1 = rfrag1.MAX() - r1pos
            }else{
                dr1 = r1pos - rfrag1.MIN()
            }

            if get_read_strand(r2) == "+"{
                dr2 = rfrag2.MAX() - r2pos
            }else{
                dr2 = r2pos - rfrag2.MIN()
            }
            fragmentsize = dr2 + dr1
        }

    return fragmentsize
}

func get_interaction_type(read1 *sam.Record, read1_chrom string, resfrag1 Interval, read2 *sam.Record,read2_chrom string, resfrag2 Interval)string{
    /*
    Returns the interaction type
    For a given reads pair and their related restriction fragment, classify
    the 3C products as :

    - Interaction
    - Self circle
    - Dangling end
    - Religation
    - Unknown
    //
    read1 = the R1 read of the pair [AlignedRead]
    read1_chrom = the chromosome of R1 read [character]
    resfrag1 = restrictin fragment overlapping the R1 read [interval]
    read2 = the R2 read of the pair [AlignedRead]
    read2_chrom = the chromosome of R2 read [character]
    resfrag2 = restrictin fragment overlapping the R2 read [interval]
    verbose = verbose mode [logical]
    */
    // If returned InteractionType=None -> Same restriction fragment
    // and same strand = Dump
    var interactionType string
 
    if !strings.Contains(r1.Flags.String(), "u") && !strings.Contains(r2.Flags.String(), "u") && resfrag1 != nil and resfrag2 != nil{
        // same restriction fragment
        if resfrag1 == resfrag2{
            // Self_circle <- ->
            if is_self_circle(read1, read2){
                interactionType = "SC"
            }else if is_dangling_end(read1, read2){
            	// Dangling_end -> <-            
                interactionType = "DE"
            }
        }else if is_religation(read1, read2, resfrag1, resfrag2){
            interactionType = "RE"
        }else{
            interactionType = "VI"
        }
    }else if !strings.Contains(r1.Flags.String(), "u") && !strings.Contains(r2.Flags.String(), "u"){
        interactionType = "SI"
    }

    return interactionType
}

/*
func get_read_tag(read, tag):
    for t in read.tags:
        if t[0] == tag:
            return t[1]
    return nil
*/

type valid struct{
    readId            string
    r1_chrom          string
    r1_pos            int
    r1_strand         int8
    r2_chrom          string
    r2_pos            int
    r2_strand         int8
    dist              int
    r1_fragname       string
    r2_fragname       string
    r1_mapq           uint8
    r2_mapq           uint8
    interactiontype   string
}

func (interacs *valid) set_valid_pair(read1, read2 *sam.Record, resFrag map[string]*Set) error{
    interacs.readId = read1.Name

    interacs.r1_chrom = read1.Ref.Name()
    interacs.r1_pos = read1.Pos
    interacs.r1_strand = get_read_strand(read1)
    interacs.r1_mapq = read1.MapQ
    interacs.r1_fragname = get_overlapping_restriction_fragment(resFrag, read1.Ref.Name(), read1)

    interacs.r2_chrom = read2.Ref.Name()
    interacs.r2_pos = read2.Pos
    interacs.r2_strand = get_read_strand(read2)
    interacs.r2_mapq = read2.MapQ
    interacs.r2_fragname = get_overlapping_restriction_fragment(resFrag, read2.Ref.Name(), read2)

    interacs.dist = get_cis_dist(read1, read2)
    return nil
}



func indeed_interaction_type(read1 *sam.Record, read2 *sam.Record, resFrag map[string]*Set, shortestInsertSize int, longestInsertSize int, minCisDist int)string{
    r1, r2 := get_ordered_reads(read1, read2)
    interacs := &valid{}
    err := interacs.set_valid_pair(r1, r2, resFrag)
    check(err)
    if interacs.r1_fragname == nil || interacs.r2_fragname == nil{
        interacs.interactiontype = "DUMP"
    }else{
        interactionType := get_interaction_type(r1, interacs.r1_chrom, interacs.r1_fragname, r2, interacs.r2_chrom, interacs.r2_fragname)
        if interactionType == "" && 



    }



}


func usage() {
    //Usage function
    fmt.Println("Usage : mapped_2hic_fragments")
    fmt.Println("-f/--fragmentFile <Restriction fragment file GFF3>")
    fmt.Println("-r/--mappedReadsFile <BAM/SAM file of mapped reads>")
    fmt.Println("[-o/--outputDir] <Output directory. Default is current directory>")
    fmt.Println("[-s/--shortestInsertSize] <Shortest insert size of mapped reads to consider>")
    fmt.Println("[-l/--longestInsertSize] <Longest insert size of mapped reads to consider>")
    fmt.Println("[-t/--shortestFragmentLength] <Shortest restriction fragment length to consider>")
    fmt.Println("[-m/--longestFragmentLength] <Longest restriction fragment length to consider>")
    fmt.Println("[-d/--minCisDist] <Minimum distance between intrachromosomal contact to consider>")
    fmt.Println("[-g/--gtag] <Genotype tag. If specified, this tag will be reported in the valid pairs output for allele specific classification>")
    fmt.Println("[-a/--all] <Write all additional output files, with information about the discarded reads (self-circle, dangling end, etc.)>")
    fmt.Println("[-S/--sam] <Output an additional SAM file with flag 'CT' for pairs classification>")
    fmt.Println("[-v/--verbose] <Verbose>")
    fmt.Println("[-h/--help] <Help>")
}

func check(e error) {
    if e != nil {
        log.Fatal(e)
    }
}

func main() {
	fragmentFile := flag.String("f", "", "fragmentFile <Restriction fragment file GFF3>")
	mappedReadsFile := flag.String("r", "", "mappedReadsFile <BAM/SAM file of mapped reads>")
	outputDir := flag.String("o", "", "<Output directory. Default is current directory>")
	shortestInsertSize := flag.Int("s", 0, "<Shortest insert size of mapped reads to consider>")
	longestInsertSize := flag.Int("l", 0, "<Longest insert size of mapped reads to consider>")
	shortestFragmentLength := flag.Int("t", 0, "<Shortest restriction fragment length to consider>")
	longestFragmentLength := flag.Int("m", 0, "<Longest restriction fragment length to consider>")
	minCisDist := flag.Int("d", 0, "<Minimum distance between intrachromosomal contact to consider>")

	_ = flag.String("g", "", "<Genotype tag. If specified, this tag will be reported in the valid pairs output for allele specific classification>")
	allOuput := flag.Bool("a", false, "<Write all additional output files, with information about the discarded reads (self-circle, dangling end, etc.)>")
	samOut := flag.Bool("S", false, "<Output an additional SAM file with flag 'CT' for pairs classification>")
	verbose := flag.Bool("v", false, "show verbose")
    help := flag.Bool("h", false, "show help")

	//stat_yes := flag.Bool("t", false, "show stat file")
	flag.Parse()

	if *help == true || *fragmentFile == "" || *mappedReadsFile == "" || *outputDir == ""{
        usage()
        log.Fatalf("Parameter error, please check -f -r -o")
    }

    fmt.Println("Genotype tag function not developed, -g parameter does not work!")

    if *verbose == true {
        fmt.Println("## overlapMapped2HiCFragments")
        fmt.Println("## mappedReadsFile=", *mappedReadsFile)
        fmt.Println("## fragmentFile=", *fragmentFile)
        fmt.Println("## minInsertSize=", *minInsertSize)
        fmt.Println("## maxInsertSize=", *maxInsertSize)
        fmt.Println("## minFragSize=", *minFragSize)
        fmt.Println("## maxFragSize=", *maxFragSize)
        fmt.Println("## allOuput=", *allOuput)
        fmt.Println("## SAM ouput=", *samOut)
        fmt.Println("## verbose=", *verbose, "\n")
    }

    // Initialize variables
    reads_counter := 0
    de_counter := 0
    re_counter := 0
    sc_counter := 0
    valid_counter := 0
    valid_counter_FF := 0
    valid_counter_RR := 0
    valid_counter_FR := 0
    valid_counter_RF := 0
    single_counter := 0
    dump_counter := 0

    /*
    // AS counter
    G1G1_ascounter := 0
    G2G2_ascounter := 0
    G1U_ascounter := 0
    UG1_ascounter := 0
    G2U_ascounter := 0
    UG2_ascounter := 0
    G1G2_ascounter := 0
    G2G1_ascounter := 0
    UU_ascounter := 0
    CF_ascounter := 0
    */

    baseReadsFile := filepath.Base(*mappedReadsFile)
    baseReadsFile = strings.TrimSuffix(strq, ".bam")
    baseReadsFile = strings.TrimSuffix(strq, ".sam")

    handle_valid, err := os.Create(fmt.Sprintf("%s/%s.validPairs", *outputDir, baseReadsFile))
    check(err)
    defer handle_valid.Close()

    if *allOuput == true {
        fmt.Println("DEPairs REPairs SCPairs ")
    }

    resFrag := load_restriction_fragment(*fragmentFile, *minFragSize, *maxFragSize, *verbose)

    bam_read, _ := os.Open(*mappedReadsFile)
    defer bam_read.Close()
    br, err := bam.NewReader(bam_read, 1)
    check(err)
    defer br.Close()

    for {
        r1, err := br.Read()
        if err == io.EOF {
            break
        }
        check(err)

        r2, err := br.Read()
        if err == io.EOF {
            break
        }
        check(err)

        if r1.Name != r2.Name {
            log.Fatalf("Forward and reverse reads not paired. Check that BAM files are sorted.")
        }

        _ = pairs(r1, r2, stat, *mapq, bw, *verbose)
    }


	bam_read1, _ := os.Open(*forward_file)
	defer bam_read1.Close()
	br1, err := bam.NewReader(bam_read1, 1)
	if err != nil {
		log.Fatalf("failed to open BAM: %v", err)
	}
	defer br1.Close()

	bam_read2, _ := os.Open(*reverse_file)
	defer bam_read2.Close()
	br2, err := bam.NewReader(bam_read2, 1)
	if err != nil {
		log.Fatalf("failed to open BAM: %v", err)
	}
	defer br2.Close()

	f_out, _ := os.OpenFile(*outfile, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0777)
	defer f_out.Close()
	bw, err := bam.NewWriter(f_out, br1.Header(), 1)
	if err != nil {
		log.Fatalf("failed to write BAM: %v", err)
	}
	defer bw.Close()

	var stat = map[string]int{
		"Total_pairs_processed":     0,
		"Unmapped_pairs":            0,
		"Pairs_with_Singleton":      0,
		"Low_qual_pairs":            0,
		"Unique_paired_alignments":  0,
		"Multiple_pairs_alignments": 0,
	}

	for {
		r1, err := br1.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("failed to read BAM1 record: %v", err)
		}

		r2, err := br2.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("failed to read BAM2 record: %v", err)
		}

		if r1.Name != r2.Name {
			log.Fatalf("Forward and reverse reads not paired. Check that BAM files are sorted.")
		}

		_ = pairs(r1, r2, stat, *mapq, bw, *verbose)
	}

	if *stat_yes {
		_ = report_stat(stat, *outfile)
	}
}
