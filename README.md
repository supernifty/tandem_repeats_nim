
## Debug
```
nim c -r main.nim < fasta_file
```

## Release
```
nim c -d:release main.nim
./main --repeat=2 --min=6 < fasta_file > bed_file
```

## Usage

* repeat: length of repeat to extract
* min: min length of repeat
