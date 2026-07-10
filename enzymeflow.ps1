param([Parameter(ValueFromRemainingArguments=$true)][string[]]$Args)
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
Push-Location $Root
try { python -m enzymeflow.cli @Args } finally { Pop-Location }

