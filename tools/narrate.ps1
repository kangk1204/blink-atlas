param([Parameter(Mandatory=$true)][string]$OutputDirectory)
$ErrorActionPreference='Stop'
$directory=(New-Item -ItemType Directory -Path $OutputDirectory -Force).FullName
$scenes=Get-Content -Raw -LiteralPath (Join-Path $PSScriptRoot 'tour.json') | ConvertFrom-Json
$voice=New-Object -ComObject SAPI.SpVoice
$english=$voice.GetVoices('Language=409')
if ($english.Count -lt 1) {throw 'An English Windows speech voice is required.'}
$voice.Voice=$english.Item(0)
$voice.Rate=-1
for ($i=0;$i -lt $scenes.Count;$i++){
  for ($j=0;$j -lt $scenes[$i].cues.Count;$j++){
    $stream=New-Object -ComObject SAPI.SpFileStream
    $stream.Open((Join-Path $directory ('{0:00}-{1:00}.wav' -f $i,$j)),3,$false)
    try {
      $voice.AudioOutputStream=$stream
      $speech=$scenes[$i].cues[$j].Replace('JF630b','J F six thirty B').Replace('JF635b','J F six thirty five B').Replace('JF629b','J F six twenty nine B').Replace('HaloTag','Halo Tag')
      [void]$voice.Speak($speech,16)
    } finally {$stream.Close()}
  }
}
