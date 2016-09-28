for y in `ls -d */`
do

  cd $y

  for x in `ls *.stderr`
  do

    grep -q "xrdcp: No such file or directory" $x

    if [ $? -eq 0 ]; then
      echo "bad job $y$x with no file error"
    fi

    grep -q "fatal" $x

    if [ $? -eq 0 ]; then
      echo "bad job $y$x with fatal error"
    fi

    grep -q "segmentation" $x

    if [ $? -eq 0 ]; then
      echo "bad job $y$x with segmentation error"
    fi

  done

  cd ./..

done

