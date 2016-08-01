sed "1s/.*/$(head -1 $1 | sed 's/"space"/"space"\tchr/g')/" $1 > temp_X221
mv temp_X221 $1;
