package main

import (
	"fmt"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"io"
	"log"
	"os"
)

func pairs(r1 *sam.Record, r2 *sam.Record, stat *map[string]int, mapq int, bw *Writer) {
	stat["Total_pairs_processed"] += 1

	if sam.IsValidRecord(r1) == false && sam.IsValidRecord(r1) == false {
		stat["Unmapped_pairs"] += 1
		return
	}

	if sam.IsValidRecord(r1) == false || sam.IsValidRecord(r1) == false {
		stat["Pairs_with_Singleton"] += 1
		return
	}

	if r1.MapQ < mapq || r2.MapQ < mapq {
		stat["Low_qual_pairs"] += 1
		return
	}

}

func main() {
	bam_read1, _ := os.Open("D:\\seqyuan\\go_bowtiePairing\\SC_Pur_01_Lib1_lane1_R1_mm9.bwt2merged.bam")
	defer bam_read1.Close()
	br1, err := bam.NewReader(bam_read1, 1)
	if err != nil {
		log.Fatalf("failed to open BAM: %v", err)
	}
	defer br1.Close()

	bam_read2, _ := os.Open("D:\\seqyuan\\go_bowtiePairing\\SC_Pur_01_Lib1_lane1_R2_mm9.bwt2merged.bam")
	defer bam_read2.Close()
	br2, err := bam.NewReader(bam_read2, 1)
	if err != nil {
		log.Fatalf("failed to open BAM: %v", err)
	}
	defer br2.Close()

	f_out, _ := os.OpenFile("D:\\seqyuan\\go_bowtiePairing\\out.bam", os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0777)
	defer f_out.Close()
	bw, err := bam.NewWriter(f_out, br.Header(), 1)
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

	aa := 0
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

		fmt.Println(r.Name, r.Ref, r.Pos, r.Seq)
		err = bw.Write(r)
		aa += 1
		if aa > 5 {
			break
		}
	}
}
