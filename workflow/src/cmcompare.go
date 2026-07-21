package main

import (
    "bufio"
    "bytes"
    "encoding/csv"
    "errors"
    "flag"
    "fmt"
    "math"
    "os"
    "os/exec"
    "path/filepath"
    "sort"
    "strconv"
    "strings"
)

type sequenceRecord struct {
    ID       string
    Sequence string
    Origin   int
    OriginName string
}

type modelSpec struct {
    Path     string
    Name     string
    Index    int
    Submodel string
    CMPath   string
}

type similarityResult struct {
    ModelA        string
    ModelB        string
    Similarity    float64
    MaximinScore  float64
    MaximinSeq    string
    Baseline      float64
}

func main() {
    var (
        k             int
        minSimilarity float64
        topK          int
        tmpdir        string
        threads       int
        evalue        float64
        outputPath    string
    )

    flag.IntVar(&k, "k", 5, "Number of sampled sequences to draw from each covariance model")
    flag.Float64Var(&minSimilarity, "min-similarity", 0.0, "Only emit pairwise similarities at or above this threshold")
    flag.IntVar(&topK, "top-k", 0, "Limit output to the top-k nearest neighbours per model (0 = all)")
    flag.StringVar(&tmpdir, "tmpdir", "", "Directory for temporary files")
    flag.IntVar(&threads, "threads", 4, "Number of threads to pass to cmsearch")
    flag.Float64Var(&evalue, "evalue", 1000.0, "Maximum E-value threshold passed to cmsearch")
    flag.StringVar(&outputPath, "output", "", "Optional TSV path for results; stdout otherwise")
    flag.Parse()

    cmPaths := flag.Args()
    if len(cmPaths) < 2 {
        fmt.Fprintf(os.Stderr, "usage: %s [flags] <cm1> <cm2> [cm3 ...]\n", filepath.Base(os.Args[0]))
        os.Exit(2)
    }

    var models []modelSpec

    if _, err := exec.LookPath("cmemit"); err != nil {
        fmt.Fprintln(os.Stderr, "error: cmemit was not found in PATH")
        os.Exit(1)
    }
    if _, err := exec.LookPath("cmsearch"); err != nil {
        fmt.Fprintln(os.Stderr, "error: cmsearch was not found in PATH")
        os.Exit(1)
    }

    workDir, err := prepareWorkDir(tmpdir)
    if err != nil {
        fmt.Fprintf(os.Stderr, "error: preparing temporary directory: %v\n", err)
        os.Exit(1)
    }
    if tmpdir == "" {
        fmt.Fprintf(os.Stderr, "temporary directory preserved at: %s\n", workDir)
    } else {
        fmt.Fprintf(os.Stderr, "using temporary directory: %s\n", workDir)
    }

    models, err = expandModels(cmPaths, workDir)
    if err != nil {
        fmt.Fprintf(os.Stderr, "error: expanding models: %v\n", err)
        os.Exit(1)
    }

    pool := make([]sequenceRecord, 0)
    seen := make(map[string]struct{})

    for i := range models {
        model := &models[i]
        resolvedCMPath, err := materializeCM(model.Path, model.Submodel, workDir, model.Index)
        if err != nil {
            fmt.Fprintf(os.Stderr, "error: materializing CM for %s: %v\n", model.Name, err)
            os.Exit(1)
        }
        model.CMPath = resolvedCMPath

        consensus, err := emitSequences(model.CMPath, true, 1)
        if err != nil {
            fmt.Fprintf(os.Stderr, "error: generating consensus for %s: %v\n", model.Name, err)
            os.Exit(1)
        }
        samples, err := emitSequences(model.CMPath, false, k)
        if err != nil {
            fmt.Fprintf(os.Stderr, "error: generating samples for %s: %v\n", model.Name, err)
            os.Exit(1)
        }

        for _, seq := range consensus {
            if _, ok := seen[seq.Sequence]; ok {
                continue
            }
            seq.ID = fmt.Sprintf("%s_consensus", sanitizeFastaID(model.Name))
            seq.Origin = model.Index
            seq.OriginName = model.Name
            pool = append(pool, seq)
            seen[seq.Sequence] = struct{}{}
        }
        for idx, seq := range samples {
            if _, ok := seen[seq.Sequence]; ok {
                continue
            }
            seq.ID = fmt.Sprintf("%s_sample_%d", sanitizeFastaID(model.Name), idx)
            seq.Origin = model.Index
            seq.OriginName = model.Name
            pool = append(pool, seq)
            seen[seq.Sequence] = struct{}{}
        }
    }

    poolPath := filepath.Join(workDir, "pool.fa")
    if err := writeFasta(poolPath, pool); err != nil {
        fmt.Fprintf(os.Stderr, "error: writing pooled FASTA: %v\n", err)
        os.Exit(1)
    }

    scoreMatrix := make(map[string][]float64, len(pool))
    for _, seq := range pool {
        scoreMatrix[seq.ID] = make([]float64, len(models))
    }

    for _, model := range models {
        tblPath := filepath.Join(workDir, fmt.Sprintf("scores_%d.tbl", model.Index))
        if err := runCmsearch(model.CMPath, poolPath, tblPath, threads, evalue); err != nil {
            fmt.Fprintf(os.Stderr, "error: scoring pool against %s: %v\n", model.Name, err)
            os.Exit(1)
        }

        knownSeqIDs := make(map[string]struct{}, len(pool))
        for _, seq := range pool {
            knownSeqIDs[seq.ID] = struct{}{}
        }

        scores, err := parseCmsearchTbl(tblPath, knownSeqIDs)
        if err != nil {
            fmt.Fprintf(os.Stderr, "error: parsing cmsearch table for %s: %v\n", model.Path, err)
            os.Exit(1)
        }

        for seqID, score := range scores {
            if row, ok := scoreMatrix[seqID]; ok {
                row[model.Index] = score
            }
        }
    }

    selfScores := make([]float64, len(models))
    for _, model := range models {
        var total float64
        count := 0
        for _, seq := range pool {
            if seq.Origin != model.Index {
                continue
            }
            total += scoreMatrix[seq.ID][model.Index]
            count++
        }
        if count > 0 {
            selfScores[model.Index] = total / float64(count)
        }
    }

    results := make([]similarityResult, 0)
    for i := 0; i < len(models); i++ {
        for j := i + 1; j < len(models); j++ {
            maximinScore := -math.MaxFloat64
            var bestSeq string
            for _, seq := range pool {
                if seq.Origin != models[i].Index && seq.Origin != models[j].Index {
                    continue
                }
                scoreA := scoreMatrix[seq.ID][i]
                scoreB := scoreMatrix[seq.ID][j]
                minScore := math.Min(scoreA, scoreB)
                if minScore > maximinScore {
                    maximinScore = minScore
                    bestSeq = seq.ID
                }
            }
            baseline := 0.0
            if len(models) > 0 {
                baseline = (selfScores[i] + selfScores[j]) / 2.0
            }
            similarity := 0.0
            if baseline != 0.0 {
                similarity = maximinScore / baseline
            }
            if similarity >= minSimilarity {
                results = append(results, similarityResult{
                    ModelA:       models[i].Name,
                    ModelB:       models[j].Name,
                    Similarity:   similarity,
                    MaximinScore: maximinScore,
                    MaximinSeq:   bestSeq,
                    Baseline:     baseline,
                })
            }
        }
    }

    sort.Slice(results, func(i, j int) bool {
        if results[i].Similarity == results[j].Similarity {
            return results[i].ModelA < results[j].ModelA
        }
        return results[i].Similarity > results[j].Similarity
    })

    if topK > 0 && len(results) > 0 {
        perModel := make(map[string][]similarityResult)
        for _, result := range results {
            perModel[result.ModelA] = append(perModel[result.ModelA], result)
        }
        filtered := make([]similarityResult, 0, len(results))
        for _, model := range models {
            modelResults := perModel[model.Name]
            if len(modelResults) == 0 {
                continue
            }
            sort.Slice(modelResults, func(i, j int) bool {
                if modelResults[i].Similarity == modelResults[j].Similarity {
                    return modelResults[i].ModelB < modelResults[j].ModelB
                }
                return modelResults[i].Similarity > modelResults[j].Similarity
            })
            if len(modelResults) > topK {
                modelResults = modelResults[:topK]
            }
            filtered = append(filtered, modelResults...)
        }
        results = filtered
    }

    var out *os.File
    if outputPath != "" {
        out, err = os.Create(outputPath)
        if err != nil {
            fmt.Fprintf(os.Stderr, "error: creating output file: %v\n", err)
            os.Exit(1)
        }
        defer out.Close()
    } else {
        out = os.Stdout
    }

    writer := csv.NewWriter(out)
    defer writer.Flush()

    _ = writer.Write([]string{"model_a", "model_b", "similarity_ratio", "empirical_maximin_score", "empirical_maximin_seq", "baseline"})
    for _, result := range results {
        _ = writer.Write([]string{
            result.ModelA,
            result.ModelB,
            fmt.Sprintf("%.6f", result.Similarity),
            fmt.Sprintf("%.6f", result.MaximinScore),
            result.MaximinSeq,
            fmt.Sprintf("%.6f", result.Baseline),
        })
    }
}

func expandModels(paths []string, workDir string) ([]modelSpec, error) {
    models := make([]modelSpec, 0, len(paths))
    for _, cmPath := range paths {
        if _, err := os.Stat(cmPath); err != nil {
            return nil, fmt.Errorf("CM file not found: %s", cmPath)
        }

        submodels, err := listModelNames(cmPath)
        if err != nil {
            return nil, err
        }

        if len(submodels) == 0 {
            resolvedPath := filepath.Join(workDir, fmt.Sprintf("model_%d.cm", len(models)))
            if err := copyFile(cmPath, resolvedPath); err != nil {
                return nil, err
            }
            models = append(models, modelSpec{
                Path:     cmPath,
                Name:     filepath.Base(cmPath),
                Index:    len(models),
                Submodel: "",
                CMPath:   resolvedPath,
            })
            continue
        }

        for _, submodel := range submodels {
            resolvedPath, err := materializeCM(cmPath, submodel, workDir, len(models))
            if err != nil {
                return nil, err
            }
            models = append(models, modelSpec{
                Path:     cmPath,
                Name:     fmt.Sprintf("%s:%s", filepath.Base(cmPath), submodel),
                Index:    len(models),
                Submodel: submodel,
                CMPath:   resolvedPath,
            })
        }
    }

    return models, nil
}

func listModelNames(cmPath string) ([]string, error) {
    cmd := exec.Command("cmstat", cmPath)
    output, err := cmd.CombinedOutput()
    if err != nil {
        return nil, fmt.Errorf("cmstat failed for %s: %w\n%s", cmPath, err, strings.TrimSpace(string(output)))
    }

    seen := make(map[string]struct{})
    names := make([]string, 0)
    for _, line := range strings.Split(string(output), "\n") {
        line = strings.TrimSpace(line)
        if line == "" || strings.HasPrefix(line, "#") || strings.HasPrefix(line, "idx") {
            continue
        }
        fields := strings.Fields(line)
        if len(fields) < 2 {
            continue
        }
        name := fields[1]
        if name == "name" || name == "acc" || name == "idx" {
            continue
        }
        if _, ok := seen[name]; ok {
            continue
        }
        seen[name] = struct{}{}
        names = append(names, name)
    }

    return names, nil
}

func materializeCM(cmPath, submodel, workDir string, idx int) (string, error) {
    targetPath := filepath.Join(workDir, fmt.Sprintf("model_%d_%s.cm", idx, sanitizeFastaID(submodel)))
    if submodel == "" {
        if err := copyFile(cmPath, targetPath); err != nil {
            return "", err
        }
        return targetPath, nil
    }

    cmd := exec.Command("cmfetch", cmPath, submodel)
    output, err := cmd.CombinedOutput()
    if err != nil {
        return "", fmt.Errorf("cmfetch failed for %s (%s): %w\n%s", cmPath, submodel, err, strings.TrimSpace(string(output)))
    }

    if err := os.WriteFile(targetPath, output, 0o644); err != nil {
        return "", err
    }
    return targetPath, nil
}

func copyFile(src, dst string) error {
    data, err := os.ReadFile(src)
    if err != nil {
        return err
    }
    return os.WriteFile(dst, data, 0o644)
}

func sanitizeFastaID(name string) string {
    replacer := strings.NewReplacer(" ", "_", "/", "_", ":", "_", "\\", "_", "\t", "_", "\n", "_", "\r", "_")
    return replacer.Replace(strings.TrimSpace(name))
}

func prepareWorkDir(tmpdir string) (string, error) {
    if tmpdir != "" {
        if err := os.MkdirAll(tmpdir, 0o755); err != nil {
            return "", err
        }
        return tmpdir, nil
    }
    return os.MkdirTemp("", "cmcompare-*")
}

func emitSequences(cmPath string, consensus bool, count int) ([]sequenceRecord, error) {
    args := []string{"-c"}
    if !consensus {
        args = []string{"-N", strconv.Itoa(count)}
    }
    args = append(args, cmPath)

    cmd := exec.Command("cmemit", args...)
    var stdout bytes.Buffer
    var stderr bytes.Buffer
    cmd.Stdout = &stdout
    cmd.Stderr = &stderr
    if err := cmd.Run(); err != nil {
        return nil, errors.New(strings.TrimSpace(stderr.String()))
    }

    return parseFasta(stdout.String(), count)
}

func parseFasta(content string, expectedCount int) ([]sequenceRecord, error) {
    scanner := bufio.NewScanner(strings.NewReader(content))
    var records []sequenceRecord
    var current sequenceRecord
    inRecord := false

    for scanner.Scan() {
        line := strings.TrimSpace(scanner.Text())
        if line == "" {
            continue
        }
        if strings.HasPrefix(line, ">") {
            if inRecord {
                records = append(records, current)
            }
            current = sequenceRecord{ID: strings.TrimPrefix(line, ">")}
            inRecord = true
            continue
        }
        if inRecord {
            current.Sequence += line
        }
    }
    if inRecord {
        records = append(records, current)
    }

    if len(records) == 0 {
        return nil, fmt.Errorf("no sequences emitted")
    }

    if expectedCount > 0 && len(records) > expectedCount {
        // cmemit can emit one consensus sequence or many sampled sequences; keep all records.
        return records, nil
    }
    return records, nil
}

func writeFasta(path string, records []sequenceRecord) error {
    f, err := os.Create(path)
    if err != nil {
        return err
    }
    defer f.Close()

    for _, record := range records {
        if _, err := fmt.Fprintf(f, ">%s\n%s\n", record.ID, record.Sequence); err != nil {
            return err
        }
    }
    return nil
}

func runCmsearch(cmPath, poolPath, tblPath string, threads int, evalue float64) error {
    cmd := exec.Command("cmsearch", "--cpu", strconv.Itoa(threads), "-E", strconv.FormatFloat(evalue, 'f', -1, 64), "--noali", "--tblout", tblPath, cmPath, poolPath)
    var stderr bytes.Buffer
    cmd.Stderr = &stderr
    if err := cmd.Run(); err != nil {
        return errors.New(strings.TrimSpace(stderr.String()))
    }
    return nil
}

func parseCmsearchTbl(path string, knownSeqIDs map[string]struct{}) (map[string]float64, error) {
    f, err := os.Open(path)
    if err != nil {
        return nil, err
    }
    defer f.Close()

    scores := make(map[string]float64)
    scanner := bufio.NewScanner(f)
    for scanner.Scan() {
        line := strings.TrimSpace(scanner.Text())
        if line == "" || strings.HasPrefix(line, "#") {
            continue
        }
        fields := strings.Fields(line)
        if len(fields) < 15 {
            continue
        }

        seqID := ""
        if _, ok := knownSeqIDs[fields[0]]; ok {
            seqID = fields[0]
        } else if len(fields) > 2 {
            if _, ok := knownSeqIDs[fields[2]]; ok {
                seqID = fields[2]
            }
        }

        if seqID == "" {
            continue
        }

        score, err := strconv.ParseFloat(fields[14], 64)
        if err != nil {
            continue
        }
        if existing, ok := scores[seqID]; !ok || score > existing {
            scores[seqID] = score
        }
    }
    if err := scanner.Err(); err != nil {
        return nil, err
    }
    return scores, nil
}
