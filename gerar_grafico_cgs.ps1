$resultado = Get-Content -Raw -Path (Join-Path $PSScriptRoot 'resultados.txt')
[System.Threading.Thread]::CurrentThread.CurrentCulture = [System.Globalization.CultureInfo]::InvariantCulture
$blocos = [regex]::Matches(
    $resultado,
    '(?ms)--- GRADIENTE CONJUGADO QUADRADO \(CGS\) ---.*?(?=\r?\n={3,}|\z)'
)

if ($blocos.Count -eq 0) {
    throw 'Nenhum resultado do CGS foi encontrado em resultados.txt.'
}

$series = @()
foreach ($blocoIndex in 0..($blocos.Count - 1)) {
    $bloco = $blocos[$blocoIndex]
    $tolerancias = @('1e-4', '1e-8')
    $tolerancia = $tolerancias[$blocoIndex]
    if ([string]::IsNullOrEmpty($tolerancia)) {
        $tolerancia = "serie $($blocoIndex + 1)"
    }

    $pontos = @{}
    foreach ($linha in [regex]::Matches($bloco.Value, 'x\[\s*(\d+)\]\s*=\s*([-+\d.eE]+)')) {
        $pontos[[int]$linha.Groups[1].Value] = [double]$linha.Groups[2].Value
    }

    if ($pontos.Count -gt 0) {
        $series += [pscustomobject]@{
            Nome = "TOL = $tolerancia"
            Pontos = $pontos
        }
    }
}

$largura = 1000
$altura = 600
$margemEsquerda = 75
$margemDireita = 30
$margemTopo = 65
$margemBase = 70
$plotLargura = $largura - $margemEsquerda - $margemDireita
$plotAltura = $altura - $margemTopo - $margemBase
$todosValores = @($series | ForEach-Object { $_.Pontos.Values })
$minimo = [math]::Floor((($todosValores | Measure-Object -Minimum).Minimum) * 10) / 10
$maximo = [math]::Ceiling((($todosValores | Measure-Object -Maximum).Maximum) * 10) / 10
if ($minimo -eq $maximo) { $maximo = $minimo + 1 }

function Escapar-Svg([string]$texto) {
    return [System.Security.SecurityElement]::Escape($texto)
}

function X([int]$indice) {
    return $margemEsquerda + ($indice / 35.0) * $plotLargura
}

function Y([double]$valor) {
    return $margemTopo + (($maximo - $valor) / ($maximo - $minimo)) * $plotAltura
}

$cores = @('#d1495b', '#00798c', '#edae49', '#30638e')
$svg = [System.Text.StringBuilder]::new()
[void]$svg.AppendLine('<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="600" viewBox="0 0 1000 600">')
[void]$svg.AppendLine('<rect width="100%" height="100%" fill="#f7f4ed"/>')
[void]$svg.AppendLine('<text x="75" y="35" font-family="Arial" font-size="22" font-weight="bold" fill="#17202a">Solucao final do CGS por componente</text>')
[void]$svg.AppendLine('<text x="75" y="55" font-family="Arial" font-size="13" fill="#52616b">Comparacao das duas tolerancias registradas em resultados.txt</text>')

for ($marca = 0; $marca -le 5; $marca++) {
    $valor = $minimo + ($maximo - $minimo) * $marca / 5
    $y = Y $valor
    [void]$svg.AppendLine(('<line x1="{0:F1}" y1="{1:F1}" x2="{2:F1}" y2="{1:F1}" stroke="#d8d2c4" stroke-width="1"/>' -f @($margemEsquerda, $y, ($largura - $margemDireita))))
    [void]$svg.AppendLine(('<text x="{0}" y="{1:F1}" text-anchor="end" font-family="Arial" font-size="12" fill="#52616b">{2:F1}</text>' -f @(($margemEsquerda - 8), ($y + 4), $valor)))
}

[void]$svg.AppendLine(('<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="#17202a" stroke-width="2"/>' -f @($margemEsquerda, ($margemTopo + $plotAltura), ($largura - $margemDireita))))
[void]$svg.AppendLine(('<line x1="{0}" y1="{1}" x2="{0}" y2="{2}" stroke="#17202a" stroke-width="2"/>' -f @($margemEsquerda, $margemTopo, ($margemTopo + $plotAltura))))

foreach ($serieIndex in 0..($series.Count - 1)) {
    $serie = $series[$serieIndex]
    $cor = $cores[$serieIndex % $cores.Count]
    $pontosSvg = @()
    for ($i = 0; $i -lt 36; $i++) {
        if ($serie.Pontos.ContainsKey($i)) {
            $pontosSvg += ('{0:F1},{1:F1}' -f (X $i), (Y $serie.Pontos[$i]))
        }
    }
    [void]$svg.AppendLine(('<polyline points="{0}" fill="none" stroke="{1}" stroke-width="3"/>' -f @(($pontosSvg -join ' '), $cor)))
    [void]$svg.AppendLine(('<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="{4}" stroke-width="4"/>' -f @(($margemEsquerda + 20 + $serieIndex * 180), 85, ($margemEsquerda + 45 + $serieIndex * 180), 85, $cor)))
    [void]$svg.AppendLine(('<text x="{0}" y="90" font-family="Arial" font-size="13" fill="#17202a">{1}</text>' -f @(($margemEsquerda + 52 + $serieIndex * 180), (Escapar-Svg $serie.Nome))))
}

foreach ($i in 0..5) {
    $indice = $i * 7
    $x = X $indice
    [void]$svg.AppendLine(('<text x="{0:F1}" y="{1}" text-anchor="middle" font-family="Arial" font-size="12" fill="#52616b">{2}</text>' -f @($x, ($altura - 42), $indice)))
}

[void]$svg.AppendLine(('<text x="{0}" y="{1}" text-anchor="middle" font-family="Arial" font-size="14" fill="#17202a">indice da incognita (i)</text>' -f ($margemEsquerda + $plotLargura / 2), ($altura - 12)))
[void]$svg.AppendLine(('<text x="18" y="{0}" transform="rotate(-90 18,{0})" text-anchor="middle" font-family="Arial" font-size="14" fill="#17202a">valor de x[i]</text>' -f ($margemTopo + $plotAltura / 2)))
[void]$svg.AppendLine('</svg>')

$caminhoSaida = Join-Path $PSScriptRoot 'grafico_cgs.svg'
[System.IO.File]::WriteAllText($caminhoSaida, $svg.ToString(), [System.Text.Encoding]::UTF8)
Write-Host "Grafico criado em $caminhoSaida"