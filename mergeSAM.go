package main

import (
	"fmt"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"io"
	"log"
	"os"
	//"strconv"
	"strings"
)

//func pairs(r1 *sam.Record, r2 *sam.Record, stat *map[string]int, mapq int, bw *Writer) {

func pairs(r1 *sam.Record, r2 *sam.Record, stat map[string]int, mapq int) error {
	v, _ := stat["Total_pairs_processed"]
	stat["Total_pairs_processed"] = v + 1

	if strings.Contains(r1.Flags.String(), "u") && strings.Contains(r2.Flags.String(), "u") {
		v, _ = stat["Unmapped_pairs"]
		stat["Unmapped_pairs"] = v + 1
		return nil
	}

	if strings.Contains(r1.Flags.String(), "u") || strings.Contains(r2.Flags.String(), "u") {
		v, _ = stat["Pairs_with_Singleton"]
		stat["Pairs_with_Singleton"] = v + 1
		return nil
	}

	if int(r1.MapQ) < mapq || int(r2.MapQ) < mapq {
		v, _ = stat["Low_qual_pairs"]
		stat["Low_qual_pairs"] = v + 1
		return nil
	}

	if is_unique_bowtie2(r1) && is_unique_bowtie2(r2) {
		v, _ = stat["Unique_paired_alignments"]
		stat["Unique_paired_alignments"] = v + 1
		return nil
	} else {
		v, _ = stat["Multiple_pairs_alignments"]
		stat["Multiple_pairs_alignments"] = v + 1
		return nil
	}

	return nil
}

func is_unique_bowtie2(r *sam.Record) (ret bool) {
	ret = false

	var primary1 int8
	var secondary1 int8

	aux, ok := r.Tag([]byte("AS"))
	if aux != nil && ok {
		aux2, ok2 := r.Tag([]byte("XS"))
		if aux2 != nil && ok2 {
			primary := aux.Value()
			secondary := aux2.Value()
			switch primary.(type) {
			case int8:
				primary1 = primary.(int8)
			case uint8:
				primary1 = int8(primary.(uint8))
			}

			switch secondary.(type) {
			case int8:
				secondary1 = secondary.(int8)
			case uint8:
				secondary1 = int8(secondary.(uint8))
			}

			if primary1 != secondary1 {
				ret = true
			}
		} else {
			ret = true
		}
	}
	return
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

	/*
		f_out, _ := os.OpenFile("D:\\seqyuan\\go_bowtiePairing\\out.bam", os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0777)
		defer f_out.Close()
		bw, err := bam.NewWriter(f_out, br.Header(), 1)
		if err != nil {
			log.Fatalf("failed to write BAM: %v", err)
		}
		defer bw.Close()
	*/

	var stat = map[string]int{
		"Total_pairs_processed":        0,
		"Unmapped_pairs":               0,
		"Pairs_with_Singleton":         0,
		"Low_qual_pairs":               0,
		"Unique_paired_alignments":     0,
		"Multiple_pairs_alignments":    0,
		"---Multiple_pairs_alignments": 0,
	}

	//aa := 0

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

		pairs(r1, r2, stat, 0)
		/*err = bw.Write(r)
		aa += 1
		if aa == 10000000 {
			break
		}
		*/
	}
	fmt.Println("Total_pairs_processed:", stat["Total_pairs_processed"])
	fmt.Println("Unmapped_pairs:", stat["Unmapped_pairs"])
	fmt.Println("Pairs_with_Singleton:", stat["Pairs_with_Singleton"])
	fmt.Println("Low_qual_pairs:", stat["Low_qual_pairs"])
	fmt.Println("Multiple_pairs_alignments:", stat["Multiple_pairs_alignments"])
	fmt.Println("Unique_paired_alignments:", stat["Unique_paired_alignments"])
}
