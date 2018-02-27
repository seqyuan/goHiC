package main

import (
	"fmt"
	"github.com/biogo/hts/bam"
	//	"github.com/biogo/hts/sam"
	"io"
	"log"
	"os"
)

func pairs(stat map[string]int32) error {
	v, _ := stat["Total_pairs_processed"]
	stat["Total_pairs_processed"] = v + 1
	//fmt.Println(stat["Total_pairs_processed"])
	return nil
}

func main() {
	bam_read, _ := os.Open("D:\\seqyuan\\go_bowtiePairing\\SC_Pur_01_Lib1_lane1_R1_mm9.bwt2merged.bam")
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

	var stat = map[string]int32{
		"Total_pairs_processed":     0,
		"Unmapped_pairs":            0,
		"Pairs_with_Singleton":      0,
		"Low_qual_pairs":            0,
		"Unique_paired_alignments":  0,
		"Multiple_pairs_alignments": 0,
	}

	_ = pairs(stat)

	aa := 0
	m2 := make(map[string]int)

	for {
		_ = pairs(stat)
		r, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("failed to read BAM record: %v", err)
		}

		err = bw.Write(r)

		aux, ok := r.Tag([]byte("AS"))
		if aux != nil && ok {
			primary := aux.Value()

			switch b := primary.(type) {
			case int:
				//fmt.Println("int", b)
				m2["int"] = 0
			case string:
				//fmt.Println("string", b)
				m2["string"] = 0
			case int8:
				//fmt.Println("int8", b)
				m2["int8"] = 0
			case uint8:
				//fmt.Println("uint8", b)
				m2["uint8"] = 0
			case int32:
				//fmt.Println("int32", b)
				m2["int32"] = 0
			default:
				//fmt.Println("else", b)
				m2["else"] = 0
				fmt.Println(b)
			}

		}

		aa += 1
		if aa > 10 {
			break
		}
	}
	fmt.Println(aa)
	fmt.Println(m2)
	fmt.Println(stat)

	qwe := int8(100)
	erio := uint8(99)

	if int8(erio) > qwe {
		fmt.Println(">")
	} else {
		fmt.Println("<")
	}
}
