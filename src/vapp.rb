
#Option parsing
def parse(args)
  parsed = {}

  args.each do |arg|
    #matches --key=value and -key=value
    match = /^-?-(?<key>.*?)(=(?<value>.*)|)$/.match(arg)
    if match
      parsed[match[:key].to_sym] = match[:value]
    else
      parsed[:text] = "#{parsed[:text]} #{arg}".strip
    end
  end

  parsed
end

def usage()
  puts "Usage: vacpp.rb -d=DIMENSION FILENAME"
end

#in ruby ARGV only contains the arguments, contrary to other languages where the first argument is the name of the script itself
if(ARGV.length()!=2)
  usage()
  exit
end
opts = parse(ARGV)
if(! opts.key?(:d)  || ! opts.key?(:text) )
  usage()
  exit
end

$ndim=opts[:d]
if(! !!($ndim =~ /[1-3]/))
  puts "Dimension should be 1,2,3"
  exit
end
inputFile = opts[:text]

#Option parsing end

$currentLineTokens = nil

$ff = File.open(inputFile)
$lineEnum = $ff.each_line
$lineNumber = 0 

$startLineIFONED    = "{^IFONED"
$startLineIFTWOD    = "{^IFTWOD"
$startLineIFTHREED  = "{^IFTHREED"
$startLineNOONED    = "{^NOONED"
$closeLine          = "}"

$LineTypeIFONED = 1
$LineTypeNOONED = 2
$LineTypeCLOSE = 3
$LineTypeFORTRAN = 4
$LineTypeIFTWOD = 5
$LineTypeIFTHREED = 6 



def nextLine()
  begin
    currentLine = ""  
    line = $lineEnum.next.rstrip
    $lineNumber = $lineNumber +1
    while(line.length() > 0 && line[-1] == "&")
      #$currentLine += line
      #replace & by space
      currentLine+=(line[0..-2] + " ")
      line = $lineEnum.next.rstrip
      $lineNumber = $lineNumber +1
    end
    currentLine+=line
    #puts("CurrentLine #{currentLine}")
    testLine = currentLine.lstrip
    if(testLine == $startLineIFONED)
      $currentLineTokens = [$LineTypeIFONED,$startLineIFONED]
    elsif(testLine == $startLineIFTWOD)
      $currentLineTokens = [$LineTypeIFTWOD,$startLineIFTWOD]
    elsif(testLine == $startLineIFTHREED)
      $currentLineTokens = [$LineTypeIFTHREED,$startLineIFTHREED]
    elsif (testLine == $startLineNOONED)
      $currentLineTokens = [$LineTypeNOONED,$startLineNOONED]
    elsif (testLine == $closeLine)
      $currentLineTokens = [$LineTypeCLOSE,$closeLine]
    else
      $currentLineTokens = [$LineTypeFORTRAN, currentLine]
    end
  rescue StopIteration
    $currentLineTokens =  nil
    $ff.close()
  end
end


def accept(str)
  return false if($currentLineTokens.nil?)
  #puts("accept  #{str} comp to #{$currentLineTokens[1]}")
  val = $currentLineTokens[1].strip
  if(val == str)
    nextLine()
    return true
  else
    return false
  end
end



def unit()
  methodArray = [method(:ifoned), method(:nooned), method(:iftwod), method(:ifthreed)]
  methodArray.each do |mm|
    r = mm.call()
    return r if(r[0])
  end
  return singleLine()
end


def singleLine()
  return [false, nil] if($currentLineTokens.nil?)
  lineType = $currentLineTokens[0]
  if(lineType == $LineTypeFORTRAN)
    retVal = [true, $currentLineTokens[1]+"\n",$lineNumber]
    nextLine()
    return retVal 
  end
  return [false, nil,$lineNumber]  
  
end


def ifstatement(acceptStr, condition)
  if(accept(acceptStr))
    r = multiunit()
    if(r[0] && accept($closeLine))
      if(condition)
        return [true,r[1],r[2]]
      else
        return [true, "",r[2]]
      end
    end
  end
  return [false, nil,$lineNumber]
end

def ifoned()
  return ifstatement($startLineIFONED,$ndim=="1")
end


def iftwod()
  return ifstatement($startLineIFTWOD,$ndim=="2")
end

def ifthreed()
  return ifstatement($startLineIFTHREED,$ndim=="3")
end


def nooned()
  return ifstatement($startLineNOONED,$ndim!="1")
end

def multiunit()
  r = unit()
  if(r[0])
      r2 = multiunit()
      if(r2[0])
        return [true, r[1] +  r2[1],r2[2]]
      else
        return [true, r[1],r[2] ]
      end  
  else
    return [false, nil,r[2]]
  end  
end


nextLine()
r = multiunit()
puts ("Syntax error on line #{r[2]}") if(!$currentLineTokens.nil?)
puts("------------------")
puts(r[1])



