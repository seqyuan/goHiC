package main

import (
	//	"bufio"
	//"flag"
	"fmt"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"io"
	"log"
	"os"
	"strings"
	"math"
)

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
	return strand
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
        r1pos = get_read_pos(r1)
        r2pos = get_read_pos(r2)
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






func main() {
	bam_file := "D:\\seqyuan\\go_bowtiePairing\\SC_Pur_01_Lib1_lane1_mm9.bwt2pairs.bam"
	bam_read1, _ := os.Open(bam_file)
	defer bam_read1.Close()
	br1, err := bam.NewReader(bam_read1, 1)
	if err != nil {
		log.Fatalf("failed to open BAM: %v", err)
	}
	defer br1.Close()

	for {
		r1, err := br1.Read()
		if err == io.EOF {
			break
		}

		r2, err := br1.Read()
		if err == io.EOF {
			break
		}

		if err != nil {
			log.Fatalf("failed to read BAM1 record: %v", err)
		}

		a := isIntraChrom(r1, r2)
		fmt.Println(a)

	}

}
