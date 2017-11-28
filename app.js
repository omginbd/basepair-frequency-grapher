/*******************
 *
 * file: stringified version of uploaded file
 *
 * returns json representation of reads:
 * {
 *  "1-mers": [
 *    {
 *      "count": 12345,
 *      "proportion": 0.75
 *    },
 *    ...
 *  ],
 *  "2-mers": [...],
 *  ...
 * }
 *
 *******************/
const parseFile = file => {}

/*******************
 *
 * read1: json representation of file1
 * read2: json representation of file2
 * (see above for data format)
 *
 * returns multidimensional array of p-values:
 * [
 *  [.1, .2, .3, ...],
 *  [...],
 *  ...
 * ]
 *
 *******************/
const compareReads = (read1, read2) => {}

const beginProcessing = () => {
  const file1Container = document.getElementById("file1")
  const file2Container = document.getElementById("file2")
  const file1Reader = new FileReader()
  const file2Reader = new FileReader()
  file1Reader.readAsText(file1Container.files[0])
  file2Reader.readAsText(file2Container.files[0])
  awaitProcessing(file1Reader, file2Reader)
}

const awaitProcessing = (f1, f2) => {
  if (f1.readyState === 2 && f2.readyState === 2) {
    const pValueTable = compareReads(parseFile(f1.result), parseFile(f2.result))
  } else {
    setTimeout(awaitProcessing, 100, f1, f2)
  }
}

// Init function
;(function() {
  const goButton = document.getElementById("go")
  goButton.addEventListener("click", beginProcessing)
})()
