filename=$(ls -1 bin/DiscreteModel{,.exe} 2>/dev/null)
extension="${filename##*.}"

if [ -n "$extension" ]; then
  extension=".$extension"
fi

cp $filename bin/`cat generated/target_name.txt`$extension
exit 0
