filename=$(ls -1 bin/OAS{,.exe} 2>/dev/null)
extension=""

if [[ "$filename" == *"OAS.exe"* ]]; then
  extension=".exe"
fi

cp $filename bin/`cat generated/target_name.txt`$extension
exit 0
