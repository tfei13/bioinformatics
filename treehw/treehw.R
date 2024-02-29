getwd

cd /Users/kylecameron/Documents/GitHub/bioinformatics/treehw

/users/kylecameron/Applications/raxml/raxml-ng \
  --all \
  --msa primate.phy \
  --model GTR+G
  --prefix T1 \
  --tree pars{10} \
  --bs-trees 200
