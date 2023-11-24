for i in {0..7}
do
    ./sol < ./org/$i.in > ./myout/$i.out
    diff -w ./myout/$i.out ./org/$i.out || break
done