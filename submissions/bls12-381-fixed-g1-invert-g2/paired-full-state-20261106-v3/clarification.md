# Editorial clarification

In checkpoint-proof.md, the phrase "full288-word comparator/control audit" means
288 instructions of word-comparison work in the predecessor, not288 data words.
Every key in this revision and its predecessor has72 words (576 bytes). The final
complete audit is paired-comparator.md. All scores and costs use36 pairs of72
words,216 comparison instructions,96 control allowance, rounded to320 per search
iteration. No different key size or alternative claim is intended.
