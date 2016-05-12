# !#/bin/bash
#
# Example:
# > ./papermovie.sh writeup/gz_hubble_data.pdf

# Inspiration: https://github.com/brownsys/paper-movies

# Issues:
#   Had to install ImageMagick (several times) with homebrew
#   Fixed font list via http://stackoverflow.com/questions/24696433/why-font-list-is-empty-for-imagemagick
#   Original codec (libx264) did not work on current version of ffmpeg

if [ "$#" -ne 1 ] || ! [ -f "$1" ]; then
  echo "Usage: $0 filename" >&2
  exit 1
fi

filename=$1

# Set the remote and/or branches for the git repository
remote=origin
branch=master

# Make sure the branch is currently up to date
git pull $remote $branch
git checkout $branch

# Find only the commits that changed this particular file
git log -- $filename  | grep "^commit" | awk '{print $2}' > commits.txt

# Download the versioned PDFs by resetting to specific commits
rev=0
while read p; do
    echo "Commit $p"
    git reset --hard $p
    if [ -f $filename ];
    then
        revpad=`printf "%03d" $rev`
        cp $filename $revpad.pdf
        
        # Use ImageMagick to turn the paginated PDF into a single image
        montage $revpad.pdf -tile 6x4 -background white -geometry 213x275-15+4 frames/$revpad.png

        let rev+=1
    fi
done < commits.txt

# Modify images in ImageMagick again to make them suitable for movie frames
# Reverse sort in this loop so that the paper goes in chronological order
rev=0
for f in `ls -1 frames/*.png | grep -v "movie\-1.png" | sort -nr`; do
    n=`printf "%03d" $rev`
    # Rename frames for ffmpeg's format style
    convert $f -crop 1280x1320+0+0 $n.png
    let rev+=1
done

# Pick the final filename
arr=$(echo $filename | tr ";" "\n")
IFS='/' read -a fnameonly <<< "$filename"
x=${fnameonly[${#fnameonly[@]}-1]}
IFS='.' read -a stub <<< "$x"

# Create the movie
ffmpeg -r 5 -i %03d.png -vcodec mpeg4 -pix_fmt yuv420p -b:v 8000k $stub.mov

# Clean up extra files in the folder
rm -f commits.txt
rm -f *.pdf
rm -f *.png
rm -f frames/*.png

# Checkout the most recent version again
git pull $remote $branch
git checkout $branch
