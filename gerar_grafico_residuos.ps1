param(
    [string]$Metodo = 'Gradiente',
    [string]$ArquivoEntrada = 'residuos_iteracoes.txt',
    [string]$ArquivoSaida = 'grafico_residuos.svg',
    [switch]$AbrirNoNavegador
)

$conteudo = Get-Content -Raw -Path (Join-Path $PSScriptRoot $ArquivoEntrada)
[System.Threading.Thread]::CurrentThread.CurrentCulture = [System.Globalization.CultureInfo]::InvariantCulture

$blocos = [regex]::Matches(
    $conteudo,
    '(?ms)Chute inicial:\s*x0\s*=\s*([-\d.]+),\s*y0\s*=\s*([-\d.]+),\s*z0\s*=\s*([-\d.]+).*?(?=\r?\nChute inicial|\z)'
)

if ($blocos.Count -eq 0) {
    throw 'Nenhum bloco de chute inicial foi encontrado em residuos_iteracoes.txt.'
}

$series = @()
foreach ($blocoIndex in 0..($blocos.Count - 1)) {
    $bloco = $blocos[$blocoIndex]
    $x0 = [double]$bloco.Groups[1].Value
    $y0 = [double]$bloco.Groups[2].Value
    $z0 = [double]$bloco.Groups[3].Value

    $pontos = [System.Collections.Generic.List[object]]::new()
    foreach ($linha in [regex]::Matches($bloco.Value, '(?m)^\s*(\d+)\s+[-\d.eE+]+\s+[-\d.eE+]+\s+[-\d.eE+]+\s+([-\d.eE+]+)\s*$')) {
        $iter = [int]$linha.Groups[1].Value
        $g = [double]$linha.Groups[2].Value
        if ($g -le 0) { $g = 1e-18 }
        $pontos.Add([pscustomobject]@{ Iter = $iter; G = $g })
    }

    if ($pontos.Count -gt 0) {
        $nome = "({0:F2}, {1:F2}, {2:F2})" -f $x0, $y0, $z0
        $series += [pscustomobject]@{
            Nome   = $nome
            Pontos = $pontos
        }
    }
}

$largura = 1100
$margemEsquerda = 85
$margemDireita = 30
$margemTopo = 70
$plotAltura = 480
$colunas = 3
$larguraLegenda = 340
$linhasLegenda = [math]::Ceiling($series.Count / $colunas)
$margemBase = 55 + 45 + ($linhasLegenda * 20)
$altura = $margemTopo + $plotAltura + $margemBase
$plotLargura = $largura - $margemEsquerda - $margemDireita

$maxIter = ($series | ForEach-Object { $_.Pontos | ForEach-Object { $_.Iter } } | Measure-Object -Maximum).Maximum
$todosLogG = @($series | ForEach-Object { $_.Pontos | ForEach-Object { [math]::Log10($_.G) } })
$minLog = [math]::Floor((($todosLogG | Measure-Object -Minimum).Minimum))
$maxLog = [math]::Ceiling((($todosLogG | Measure-Object -Maximum).Maximum))
if ($minLog -eq $maxLog) { $maxLog = $minLog + 1 }

function X([int]$iter) {
    return $margemEsquerda + ($iter / $maxIter) * $plotLargura
}

function Y([double]$g) {
    $logG = [math]::Log10($g)
    return $margemTopo + (($maxLog - $logG) / ($maxLog - $minLog)) * $plotAltura
}

function Escapar-Svg([string]$texto) {
    return [System.Security.SecurityElement]::Escape($texto)
}

$cores = @('#d1495b', '#00798c', '#edae49', '#30638e', '#8a5a44', '#7a9d54', '#9b5de5', '#f15bb5', '#00bbf9', '#fee440', '#264653')
$svg = [System.Text.StringBuilder]::new()
[void]$svg.AppendLine(('<svg xmlns="http://www.w3.org/2000/svg" width="{0}" height="{1}" viewBox="0 0 {0} {1}">' -f $largura, $altura))
[void]$svg.AppendLine('<rect width="100%" height="100%" fill="#f7f4ed"/>')
[void]$svg.AppendLine(( '<text x="85" y="32" font-family="Arial" font-size="22" font-weight="bold" fill="#17202a">Convergencia do Metodo de {0} por chute inicial</text>' -f [System.Security.SecurityElement]::Escape($Metodo)))
[void]$svg.AppendLine(( '<text x="85" y="52" font-family="Arial" font-size="13" fill="#52616b">G(x,y,z) (escala log10) por iteracao - {0}</text>' -f [System.Security.SecurityElement]::Escape($ArquivoEntrada)))

for ($marca = 0; $marca -le ($maxLog - $minLog); $marca++) {
    $valorLog = $minLog + $marca
    $y = Y ([math]::Pow(10, $valorLog))
    [void]$svg.AppendLine(('<line x1="{0:F1}" y1="{1:F1}" x2="{2:F1}" y2="{1:F1}" stroke="#d8d2c4" stroke-width="1"/>' -f @($margemEsquerda, $y, ($largura - $margemDireita))))
    [void]$svg.AppendLine(('<text x="{0}" y="{1:F1}" text-anchor="end" font-family="Arial" font-size="12" fill="#52616b">1e{2}</text>' -f @(($margemEsquerda - 8), ($y + 4), $valorLog)))
}

[void]$svg.AppendLine(('<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="#17202a" stroke-width="2"/>' -f @($margemEsquerda, ($margemTopo + $plotAltura), ($largura - $margemDireita))))
[void]$svg.AppendLine(('<line x1="{0}" y1="{1}" x2="{0}" y2="{2}" stroke="#17202a" stroke-width="2"/>' -f @($margemEsquerda, $margemTopo, ($margemTopo + $plotAltura))))

foreach ($serieIndex in 0..($series.Count - 1)) {
    $serie = $series[$serieIndex]
    $cor = $cores[$serieIndex % $cores.Count]
    $pontosSvg = $serie.Pontos | ForEach-Object { '{0:F1},{1:F1}' -f (X $_.Iter), (Y $_.G) }
    [void]$svg.AppendLine(('<polyline points="{0}" fill="none" stroke="{1}" stroke-width="2"/>' -f @(($pontosSvg -join ' '), $cor)))
}

foreach ($serieIndex in 0..($series.Count - 1)) {
    $serie = $series[$serieIndex]
    $cor = $cores[$serieIndex % $cores.Count]
    $coluna = $serieIndex % $colunas
    $linha = [math]::Floor($serieIndex / $colunas)
    $lx = $margemEsquerda + $coluna * $larguraLegenda
    $ly = $margemTopo + $plotAltura + 45 + $linha * 20
    [void]$svg.AppendLine(('<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="{3}" stroke-width="4"/>' -f @($lx, $ly, ($lx + 25), $cor)))
    [void]$svg.AppendLine(('<text x="{0}" y="{1}" font-family="Arial" font-size="12" fill="#17202a">{2}</text>' -f @(($lx + 32), ($ly + 4), (Escapar-Svg $serie.Nome))))
}

foreach ($i in 0..5) {
    $iter = [math]::Round($maxIter * $i / 5)
    $x = X $iter
    [void]$svg.AppendLine(('<text x="{0:F1}" y="{1}" text-anchor="middle" font-family="Arial" font-size="12" fill="#52616b">{2}</text>' -f @($x, ($margemTopo + $plotAltura + 18), $iter)))
}

[void]$svg.AppendLine(('<text x="{0}" y="{1}" text-anchor="middle" font-family="Arial" font-size="14" fill="#17202a">iteracao</text>' -f ($margemEsquerda + $plotLargura / 2), ($margemTopo + $plotAltura + 34)))
[void]$svg.AppendLine(('<text x="20" y="{0}" transform="rotate(-90 20,{0})" text-anchor="middle" font-family="Arial" font-size="14" fill="#17202a">G(x,y,z) (log10)</text>' -f ($margemTopo + $plotAltura / 2)))
[void]$svg.AppendLine('</svg>')

$caminhoSaida = Join-Path $PSScriptRoot $ArquivoSaida
[System.IO.File]::WriteAllText($caminhoSaida, $svg.ToString(), [System.Text.Encoding]::UTF8)
Write-Host "Grafico criado em $caminhoSaida"

if ($AbrirNoNavegador) {
    $uriSaida = [System.Uri]::new([System.IO.Path]::GetFullPath($caminhoSaida)).AbsoluteUri
    try {
        Start-Process -FilePath $uriSaida -ErrorAction Stop
    }
    catch {
        Write-Error "Nao foi possivel abrir o grafico no navegador: $_"
        exit 1
    }
}
