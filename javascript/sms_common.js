//    Sequence Manipulation Suite. A collection of simple JavaScript programs
//    for generating, formatting, and analyzing short DNA and protein
//    sequences.
//    Copyright (C) 2020 Paul Stothard stothard@ualberta.ca
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

//Written by Paul Stothard, University of Alberta, Canada

function checkGeneticCode(arrayOfPatterns) {
  var z = 0;
  var codon = "";
  var oneMatch = false;
  var testSequence =
    "gggggaggtggcgaggaagatgacgtggtagttgtcgcggcagctgccaggagaagtagcaagaaaaataacatgataattatcacgacaactacctggtgatgttgctagtaatattacttgttatttttctcgtcatcttcccggcgacgtcgccagcaacatcacctgctacttctcccgccacctccc";
  while (z < arrayOfPatterns.length) {
    if (arrayOfPatterns[z].search(/^\s*\/[a-zA-Z\|\[\]]+\/=[a-zA-Z\*]/) == -1) {
      alert(
        "Genetic code error: one or more patterns have been entered incorrectly."
      );
      return false;
    }
    if (moreExpressionCheck(arrayOfPatterns[z]) == false) {
      alert(
        "Genetic code error: one or more patterns have been entered incorrectly."
      );
      return false;
    }
    z = z + 1;
  }
  var geneticCodeMatchResult = new Array(arrayOfPatterns.length);
  var geneticCodeMatchExp = new Array(arrayOfPatterns.length);
  for (var j = 0; j < arrayOfPatterns.length; j++) {
    geneticCodeMatchExp[j] = eval(arrayOfPatterns[j].match(/\/.+\//) + "gi");
    geneticCodeMatchResult[j] = arrayOfPatterns[j]
      .match(/=[a-zA-Z\*]/)
      .toString();
    geneticCodeMatchResult[j] = geneticCodeMatchResult[j].replace(/=/g, "");
  }
  for (var i = 0; i <= testSequence.length - 3; i = i + 3) {
    codon = testSequence.substring(i, i + 3);
    for (var j = 0; j < geneticCodeMatchExp.length; j++) {
      if (codon.search(geneticCodeMatchExp[j]) != -1) {
        if (oneMatch == true) {
          alert(
            "Genetic code error: more than one amino acid is coded by the codon: " +
              codon +
              "."
          );
          return false;
        }
        oneMatch = true;
      }
    }
    if (oneMatch == false) {
      alert("The genetic code expressions are missing a codon.");
      return false;
    }
    oneMatch = false;
  }
  return true;
}

function checkRestPatterns(arrayOfPatterns) {
  var z = 0;
  while (z < arrayOfPatterns.length) {
    if (
      arrayOfPatterns[z].search(
        /^\s*\/[acgturyswkmbdhvn\[\]]+\/\s+\([^\/]+\)\d+/i
      ) == -1
    ) {
      alert("One or more patterns have been entered incorrectly.");
      return false;
    }
    if (moreExpressionCheck(arrayOfPatterns[z]) == false) {
      alert("One or more patterns have been entered incorrectly.");
      return false;
    }
    z = z + 1;
  }
  return true;
}

function closeWindow() {
  outputWindow.document.write("</body>\n</html>\n");
  outputWindow.status = "Done.";
  outputWindow.document.close();
  return true;
}

function convertDegenerates(sequence) {
  sequence = sequence.toLowerCase();
  sequence = sequence.replace(/t/g, "[TU]");
  sequence = sequence.replace(/r/g, "[AGR]");
  sequence = sequence.replace(/y/g, "[CTUY]");
  sequence = sequence.replace(/s/g, "[GCS]");
  sequence = sequence.replace(/w/g, "[ATUW]");
  sequence = sequence.replace(/k/g, "[GTUK]");
  sequence = sequence.replace(/m/g, "[ACM]");
  sequence = sequence.replace(/b/g, "[CGTUBSKY]");
  sequence = sequence.replace(/d/g, "[AGTUDRKW]");
  sequence = sequence.replace(/h/g, "[ACTUHMYW]");
  sequence = sequence.replace(/v/g, "[ACGVSMR]");
  sequence = sequence.replace(/n/g, "[ACGTURYSWKMBDHVN]");
  return sequence;
}

function getArrayOfFasta(sequenceData) {
  var arrayOfFasta = new Array();
  var matchArray;
  var re = /\>[^\>]+/g;
  if (sequenceData.search(/\>[^\f\n\r]+[\f\n\r]/) != -1) {
    while ((matchArray = re.exec(sequenceData))) {
      arrayOfFasta.push(matchArray[0]);
    }
  } else {
    arrayOfFasta[0] = sequenceData;
  }
  return arrayOfFasta;
}

function getFastaTitleFromTitleAndSequence(fastaSequenceTitle, sequence) {
  var stringToReturn =
    "&gt;results for " + sequence.length + " residue sequence ";
  if (fastaSequenceTitle.search(/[^\s]/) != -1) {
    stringToReturn = stringToReturn + '"' + fastaSequenceTitle + '"';
  }
  stringToReturn =
    stringToReturn + ' starting "' + sequence.substring(0, 10) + '"';
  return stringToReturn + "\n";
}

function getGeneticCodeMatchExp(arrayOfPatterns) {
  var geneticCodeMatchExp = new Array(arrayOfPatterns.length);
  for (var j = 0; j < arrayOfPatterns.length; j++) {
    geneticCodeMatchExp[j] = eval(arrayOfPatterns[j].match(/\/.+\//) + "gi");
  }
  return geneticCodeMatchExp;
}

function getGeneticCodeMatchResult(arrayOfPatterns) {
  var geneticCodeMatchResult = new Array(arrayOfPatterns.length);
  for (var j = 0; j < arrayOfPatterns.length; j++) {
    geneticCodeMatchResult[j] = arrayOfPatterns[j]
      .match(/=[a-zA-Z\*]/)
      .toString();
    geneticCodeMatchResult[j] = geneticCodeMatchResult[j].replace(/=/g, "");
  }
  return geneticCodeMatchResult;
}

function getInfoFromTitleAndSequence(fastaSequenceTitle, sequence) {
  var stringToReturn = "Results for " + sequence.length + " residue sequence ";
  if (fastaSequenceTitle.search(/[^\s]/) != -1) {
    stringToReturn = stringToReturn + '"' + fastaSequenceTitle + '"';
  }
  stringToReturn =
    stringToReturn + ' starting "' + sequence.substring(0, 10) + '"';
  return '<div class="info">' + stringToReturn + "</div>\n";
}

function getSequenceFromFasta(sequenceRecord) {
  if (sequenceRecord.search(/\>[^\f\n\r]+[\f\n\r]/) != -1) {
    sequenceRecord = sequenceRecord.replace(/\>[^\f\n\r]+[\f\n\r]/, "");
  }
  return sequenceRecord;
}

function getTitleFromFasta(sequenceRecord) {
  var fastaTitle = "Untitled";
  if (sequenceRecord.search(/\>[^\f\n\r]+[\f\n\r]/) != -1) {
    fastaTitle = sequenceRecord.match(/\>[^\f\n\r]+[\f\n\r]/, "").toString();
    fastaTitle = fastaTitle.replace(/\>|[\f\n\r]/g, "");
    fastaTitle = fastaTitle.replace(/\s{2,}/g, " ");
    fastaTitle = fastaTitle.replace(/[\<\>]/gi, "");
  }
  return fastaTitle;
}

function moreExpressionCheck(expressionToCheck) {
  if (
    expressionToCheck.search(/\[[A-Za-z\|]*\[/) != -1 ||
    expressionToCheck.search(/\][A-Za-z\|]*\]/) != -1 ||
    expressionToCheck.search(/\[\]/) != -1 ||
    expressionToCheck.search(/\/[A-Za-z\|]*\]/) != -1 ||
    expressionToCheck.search(/\[[A-Za-z\|]*\//) != -1 ||
    expressionToCheck.search(/\|\|/) != -1 ||
    expressionToCheck.search(/\/\|/) != -1 ||
    expressionToCheck.search(/\|\//) != -1 ||
    expressionToCheck.search(/\[.\]/) != -1 ||
    expressionToCheck.search(/\</) != -1 ||
    expressionToCheck.search(/\>/) != -1
  ) {
    return false;
  }
  return true;
}

function openWindow(title) {
  _openWindow(title, true);
}

function _openWindow(title, isColor) {
  outputWindow = window.open(
    "",
    "my_new_window",
    "toolbar=no, location=no, directories=no, status=yes, menubar=yes, scrollbars=yes, resizable=yes, copyhistory=no, width=800, height=400"
  );
  outputWindow.focus();
  outputWindow.document.write(
    '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n' +
      '<html lang="en">\n' +
      "<head>\n" +
      "<title>Sequence Manipulation Suite</title>\n" +
      '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1" />\n'
  );

  if (isColor) {
    outputWindow.document.write(
      '<style type="text/css">\n' +
        "body.main {font-size:90%; font-family: arial, sans-serif; color: #000000; background-color: #FFFFFF}\n" +
        "div.pre {color: #000000; font-family: courier, sans-serif; white-space: pre}\n" +
        "div.title {color: #000000; text-align: left; background-color: #FFFFFF}\n" +
        "div.info {font-weight: bold}\n" +
        "span.none, td.none {color: #000000; background-color: #FFFFFF}\n" +
        "span.one, td.one {color: #000000; background-color: #66FF00}\n" +
        "span.two, td.two {color: #000000; background-color: #FFFF66}\n" +
        "span.three, td.three {color: #000000; background-color: #FFFFFF}\n" +
        "span.forward_primer, td.forward_primer {color: #000000; background-color: #FF66FF}\n" +
        "span.reverse_primer, td.reverse_primer {color: #000000; background-color: #FF9933}\n" +
        "span.current_sequence {color: #000000; background-color: #FFFFFF}\n" +
        "span.mutated_sequence {color: #990066; background-color: #FFFFFF}\n" +
        "td.many {color: #000000}\n" +
        "td.title {font-weight: bold; color: #000000; background-color: #FFFFFF}\n" +
        "</style>\n"
    );
  } else {
    outputWindow.document.write(
      '<style type="text/css">\n' +
        "body.main {font-size:90%; font-family: arial, sans-serif; color: #000000; background-color: #FFFFFF; margin: 0 auto; padding: 0}\n" +
        "div.pre {color: #000000; background-color: #FFFFFF; font-family: courier, sans-serif; white-space: pre}\n" +
        "div.title {display: none}\n" +
        "div.info {font-weight: bold}\n" +
        "span.none, td.none {color: #000000; background-color: #FFFFFF}\n" +
        "span.one, td.one {color: #000000; text-decoration: underline; background-color: #FFFFFF}\n" +
        "span.two, td.two {color: #000000; font-style: italic; background-color: #FFFFFF}\n" +
        "span.three, td.three {color: #000000; background-color: #FFFFFF}\n" +
        "span.forward_primer, td.forward_primer {color: #000000; background-color: #FFFFFF}\n" +
        "span.reverse_primer, td.reverse_primer {color: #000000; background-color: #FFFFFF}\n" +
        "span.current_sequence {color: #000000; background-color: #FFFFFF}\n" +
        "span.mutated_sequence {color: #000000; text-decoration: underline; background-color: #FFFFFF}\n" +
        "td.many {color: #000000; background-color: #FFFFFF}\n" +
        "td.title {font-weight: bold; color: #000000; background-color: #FFFFFF}\n" +
        "img {display: none}\n" +
        "</style>\n"
    );
  }
  outputWindow.document.write(
    "</head>\n" +
      '<body class="main">\n' +
      '<div class="title">' +
      title +
      " results</div>\n"
  );
  outputWindow.status = "Please Wait.";
  return true;
}

function openWindowAlign(title) {
  _openWindowAlign(title, true);
}

function _openWindowAlign(title, isBackground) {
  outputWindow = window.open(
    "",
    "my_new_window",
    "toolbar=no, location=no, directories=no, status=yes, menubar=yes, scrollbars=yes, resizable=yes, copyhistory=no, width=800, height=400"
  );
  outputWindow.focus();
  outputWindow.document.write(
    '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n' +
      '<html lang="en">\n' +
      "<head>\n" +
      "<title>Sequence Manipulation Suite</title>\n" +
      '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1" />\n'
  );
  if (isBackground) {
    outputWindow.document.write(
      '<style type="text/css">\n' +
        "body.main {font-family: arial, sans-serif; font-size:90%; color: #000000; background-color: #FFFFFF}\n" +
        "div.pre {color: #000000; font-family: courier, sans-serif; white-space: pre}\n" +
        "div.title {color: #000000; text-align: left; background-color: #FFFFFF}\n" +
        "div.info {font-weight: bold}\n" +
        "span.ident {color: #FFFFFF; background-color: #000000}\n" +
        "span.sim {color: #FFFFFF; background-color: #666666}\n" +
        "span.g, span.a, span.v, span.l, span.i {color: #000000; background-color: #C0C0C0}\n" +
        "span.f, span.y, span.w {color: #000000; background-color: #FF6600}\n" +
        "span.c, span.m {color: #000000; background-color: #FFFF00}\n" +
        "span.s, span.t {color: #000000; background-color: #66FF00}\n" +
        "span.k, span.r, span.h {color: #000000; background-color: #FF0000}\n" +
        "span.d, span.e {color: #000000; background-color: #0066FF}\n" +
        "span.n, span.q {color: #000000; background-color: #996633}\n" +
        "span.p {color: #000000; background-color: #FF99FF}\n" +
        "</style>\n"
    );
  } else {
    outputWindow.document.write(
      '<style type="text/css">\n' +
        "body.main {font-family: arial, sans-serif; font-size:90%; color: #000000; background-color: #FFFFFF}\n" +
        "div.pre {color: #000000; font-family: courier, sans-serif; white-space: pre}\n" +
        "div.title {display: none}\n" +
        "div.info {font-weight: bold}\n" +
        "span.ident {color: #000000; font-weight: bold; text-decoration: underline; background-color: #FFFFFF}\n" +
        "span.sim {color: #000000; font-weight: bold; background-color: #FFFFFF}\n" +
        "span.diff {color: #999999; background-color: #FFFFFF}\n" +
        "span.g, span.a, span.v, span.l, span.i {color: #CC33CC; background-color: #FFFFFF}\n" +
        "span.f, span.y, span.w {color: #FF6600; background-color: #FFFFFF}\n" +
        "span.c, span.m {color: #FFCC00; background-color: #FFFFFF}\n" +
        "span.s, span.t {color: #CCFF00; background-color: #FFFFFF}\n" +
        "span.k, span.r, span.h {color: #FF0000; background-color: #FFFFFF}\n" +
        "span.d, span.e {color: #0000FF; background-color: #FFFFFF}\n" +
        "span.n, span.q {color: #996633; background-color: #FFFFFF}\n" +
        "span.p {color: #00FFCC; background-color: #FFFFFF}\n" +
        "img {display: none}\n" +
        "</style>\n"
    );
  }
  outputWindow.document.write(
    "</head>\n" +
      '<body class="main">\n' +
      '<div class="title">' +
      title +
      " results</div>\n"
  );
  outputWindow.status = "Please Wait.";
  return true;
}

function removeNonDna(sequence) {
  return sequence.replace(/[^gatucryswkmbdhvnxGATUCRYSWKMBDHVNX]/g, "");
}

function rightNum(theNumber, sequenceToAppend, lengthOfColumn, tabIn) {
  var j = 0;
  var tempString = "";
  theNumber = theNumber.toString();
  for (var j = theNumber.length; j < lengthOfColumn; j++) {
    tempString = tempString + " ";
  }
  theNumber = tempString + theNumber + " ";
  sequenceToAppend = sequenceToAppend + theNumber + tabIn;
  return sequenceToAppend;
}

function testScript() {
  //test some javascript functions to see how the browser performs.
  //want to prevent non Javascript 1.5 browsers from attempting to run.
  //first test Array.push()
  var testArray = new Array();
  var testString = "1234567890";

  testArray.push(testString);
  if (testArray[0] != testString) {
    alert(
      "Array object push method not supported. See browser compatibility page."
    );
    return false;
  }

  //now test the 'm' flag in a regular expression
  testString = "1\n2\n3";
  var re = /^2$/m;
  if (testString.search(re) == -1) {
    alert(
      "Regular expression 'm' flag not supported. See browser compatibility page."
    );
    return false;
  }

  var caughtException = false;
  //now test exception handling
  try {
    re = eval(
      "Exception handling not supported. Check browser compatibility page."
    );
  } catch (e) {
    caughtException = true;
  }

  if (!caughtException) {
    alert("Exception handling not supported. See browser compatibility page.");
  }

  //now test replace lambda function
  testString = "123";
  testString = testString.replace(/(\d)/g, function (str, p1, offset, s) {
    return p1 + "X";
  });
  if (testString != "1X2X3X") {
    alert(
      "Nested function in String replace method not supported. See browser compatibility page."
    );
    return false;
  }

  //test number methods toFixed() and toPrecision()
  var testNum = 2489.8237;
  if (testNum.toFixed(3) != 2489.824) {
    alert(
      "Number toFixed() method not supported. See browser compatibility page."
    );
    return false;
  }

  if (testNum.toPrecision(5) != 2489.8) {
    alert(
      "Number toPrecision() method not supported. See browser compatibility page."
    );
    return false;
  }

  return true;
}

function writeRestrictionSites(sequence, arrayOfItems, dnaConformation) {
  var resultArray = new Array();
  var lookAhead = 50;
  var lowerLimit = 0;
  var upperLimit = sequence.length;
  var shiftValue = 0;
  var cutDistance;
  var matchExp;
  var matchPosition;
  var tempString;
  var backGroundClass;
  var matchArray;
  var timesFound = 0;
  if (dnaConformation == "circular") {
    shiftValue = sequence.substring(0, lookAhead).length;
    sequence =
      sequence.substring(sequence.length - lookAhead, sequence.length) +
      sequence +
      sequence.substring(0, lookAhead);
    lowerLimit = 0 + shiftValue;
    upperLimit = upperLimit + shiftValue;
  }
  outputWindow.document.write(
    '<table border="1" width="100%" cellspacing="0" cellpadding="2"><tbody>\n'
  );
  outputWindow.document.write(
    '<tr><td class="title" width="200px">' +
      "Site:" +
      '</td><td class="title">' +
      "Positions:" +
      "</td></tr>\n"
  );
  for (var i = 0; i < arrayOfItems.length; i++) {
    tempString = "none";
    backGroundClass = "many";
    matchExp = arrayOfItems[i].match(/\/.+\//) + "gi";
    matchPosition = 0;
    matchExp = eval(matchExp);
    cutDistance = parseFloat(
      arrayOfItems[i]
        .match(/\)\D*\d+/)
        .toString()
        .replace(/\)\D*/, "")
    );

    while ((matchArray = matchExp.exec(sequence))) {
      matchPosition = matchExp.lastIndex - cutDistance;
      if (matchPosition >= lowerLimit && matchPosition < upperLimit) {
        timesFound++;
        tempString = tempString + ", " + (matchPosition - shiftValue + 1);
      }
      matchExp.lastIndex = matchExp.lastIndex - RegExp.lastMatch.length + 1;
    }

    if (tempString.search(/\d/) != -1) {
      tempString = tempString.replace(/none,\s*/, "");
    }

    if (timesFound == 0) {
      backGroundClass = "none";
    } else if (timesFound == 1) {
      backGroundClass = "one";
    } else if (timesFound == 2) {
      backGroundClass = "two";
    } else if (timesFound == 3) {
      backGroundClass = "three";
    } else {
      backGroundClass = "many";
    }

    outputWindow.document.write(
      '<tr><td class="' +
        backGroundClass +
        '">' +
        arrayOfItems[i]
          .match(/\([^\(]+\)/)
          .toString()
          .replace(/\(|\)/g, "") +
        '</td><td class="' +
        backGroundClass +
        '">' +
        tempString +
        "</td></tr>\n"
    );

    timesFound = 0;
  }
  outputWindow.document.write("</tbody></table>\n");
  return true;
}
