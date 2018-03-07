package main

import (
	//	"bufio"
	"flag"
	"fmt"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"io"
	"log"
	"os"
	"strings"
)

func usage() {
	fmt.Println("Usage : ./mergeSAM")
	fmt.Println("-f <forward read mapped file>")
	fmt.Println("-r <reverse read mapped file>")
	fmt.Println("-o <Output file>")
	fmt.Println("[-q] <minimum reads mapping quality,default 0>")
	fmt.Println("[-t] <generate a stat file>")
	fmt.Println("[-v] <Verbose>")
}

func pairs(r1 *sam.Record, r2 *sam.Record, stat map[string]int, mapq int, bw *bam.Writer, verbose bool) error {
	v, _ := stat["Total_pairs_processed"]
	stat["Total_pairs_processed"] = v + 1
	if verbose && v%1000000 == 0 {
		fmt.Println("##", v)
	}

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
		err := bw.Write(r1)
		if err != nil {
			log.Fatalf("write outbam err")
		}
		err = bw.Write(r2)
		if err != nil {
			log.Fatalf("write outbam err")
		}
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

func percent(part int, all int) string {
	percent := fmt.Sprintf("%.3f", float64(part)*100/float64(all))
	return percent
}

func report_stat(stat map[string]int, outfile string) error {
	out_stat := strings.Replace(outfile, ".bam", ".pairstat", -1)
	//stat_out, _ := os.OpenFile(out_stat, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0777)
	stat_out, _ := os.Create(out_stat)
	defer stat_out.Close()
	//writer := bufio.NewWriter(stat_out)

	keys := [6]string{"Total_pairs_processed", "Unmapped_pairs", "Pairs_with_Singleton", "Low_qual_pairs", "Unique_paired_alignments", "Multiple_pairs_alignments"}
	for _, k := range keys {
		outbyte := fmt.Sprintf("%v\t%v\t%v\n", k, stat[k], percent(stat[k], stat["Total_pairs_processed"]))
		_, _ = stat_out.WriteString(outbyte)
	}
	//writer.Flush()
	return nil
}

//%v\t%v  k, stat[k],

func main() {
	forward_file := flag.String("f", "", "forward read mapped file")
	reverse_file := flag.String("r", "", "reverse read mapped file")
	outfile := flag.String("o", "", "out put file")
	mapq := flag.Int("q", 0, "minimum reads mapping quality")
	stat_yes := flag.Bool("t", false, "show stat file")
	verbose := flag.Bool("v", false, "show verbose")
	flag.Parse()

	if *forward_file == "" || *reverse_file == "" || *outfile == "" {
                usage()
		log.Fatalf("Parameter error,please check -f -r -o")
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

