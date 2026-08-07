# PowerShell script to copy parser sources into the project.
# Working directory should be project root.
# This does not currently update the main library sources, which would come from
# https://github.com/tree-sitter/tree-sitter.

if ((Get-Item .).Name -ne 'turbowave') {
    Write-Error 'Must be run from TW project root'
    exit(1)
}

Add-Type -AssemblyName System.Windows.Forms
$selector = New-Object System.Windows.Forms.FolderBrowserDialog
$selector.Description = "Select TS project"

if ($selector.ShowDialog() -eq "OK") {
    $src = $selector.SelectedPath
    $cdest = "core\source\io\parser.c"
    $hdest = "core\source\io\tree-sitter\lib\include\tree_sitter\parser.h"
    Copy-Item -Path $src\src\parser.c -Destination $cdest
    Copy-Item -Path $src\src\tree_sitter\parser.h -Destination $hdest
}
