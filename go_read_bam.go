package main

import (
	"fmt"
	"github.com/biogo/hts/bam"
	"io"
	"log"
	"os"
)

func main() {
	bam_file, _ := os.Open("D:\\gobismark\\go_read_bam_test\\b1_Lib1_lane1_hg19_pure.bwt2pairs.bam")
	// Create a BAI for the BAM read from standard in and write it to standard out.
	br, err := bam.NewReader(bam_file, 1)
	if err != nil {
		log.Fatalf("failed to open BAM: %v", err)
	}
	aa := 0
	for {
		r, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("failed to read BAM record: %v", err)
		}

		fmt.Println(r.Name, r.Ref, r.Pos, r.Seq)
		aa += 1
		if aa > 3 {
			break
		}
	}
}
