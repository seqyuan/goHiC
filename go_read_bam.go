package main

import (
	"fmt"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"io"
	"log"
	"os"
)

func main() {
	bam_read, _ := os.Open("D:\\seqyuan\\go_bowtiePairing\\SC_Pur_01_Lib1_lane1_mm9.bwt2pairs.bam")
	defer bam_read.Close()
	br, err := bam.NewReader(bam_read, 1)
	if err != nil {
		log.Fatalf("failed to open BAM: %v", err)
	}
	defer br.Close()

	f_out, _ := os.OpenFile("D:\\seqyuan\\go_bowtiePairing\\out.bam", os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0777)
	defer f_out.Close()
	bw, err := bam.NewWriter(f_out, br.Header(), 1)
	if err != nil {
		log.Fatalf("failed to write BAM: %v", err)
	}
	defer bw.Close()

	aa := 0
	for {
		r, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("failed to read BAM record: %v", err)
		}

		_, ok := r.Tag([]byte("ASsss"))

		fmt.Println(r.Name, ok)
		if sam.IsValidRecord(r) == false {
			fmt.Println("000000000000000")
		}

		err = bw.Write(r)
		aa += 1
		if aa > 5 {
			break
		}
	}
}
