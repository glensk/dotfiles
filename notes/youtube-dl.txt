you have to quote the url! (at least in zsh)
youtube-dl -x --audio-format mp3 'https://www.youtube.com/watch?v=58KS5nsnM2g'
    -x, --extract-audio              Convert video files to audio-only files (requires ffmpeg or avconv and ffprobe or avprobe)

for videoad use -f 22 to get best quality and without " "
youtube-dl -f 22 https://www.youtube.com/watch\?v\=NAmBqEBI2NI
