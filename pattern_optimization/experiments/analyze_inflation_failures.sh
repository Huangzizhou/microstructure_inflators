# Parses through pattern optimization output(s) and summarizes the causes of
# inflation failure.
grep -h Except $@ -A2 | awk "NR%4 == 3" | sort | uniq -c
