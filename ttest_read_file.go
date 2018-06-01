package main

import (
	//	"bufio"
	"bufio"
	"fmt"
	"io"
	//"log"
	"os"
)

func read1(file string) {
	rw, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer rw.Close()
	rb := bufio.NewReader(rw)
	for {
		line, _, err := rb.ReadLine()
		if err == io.EOF {
			break
		}
		//do something
		for _, l := range bytes.Split(text, []byte{'\t'}) {
			fmt.Println(string(l[5]))
		}
	}
}

func read2(file string) {
	rw, err := os.Open("")
	if err != nil {
		panic(err)
	}
	defer rw.Close()
	sb := bufio.NewScanner(rw)
	for sb.Scan() {
		//do something
		fmt.Println(sb.Text())
	}
	if err := sb.Err(); err != nil {
		panic(err)
	}
}

func main() {
	userFile := "D:\\TMK\\HindIII_resfrag_Arabidopsis_thaliana.bed"

	read1(userFile)
}
