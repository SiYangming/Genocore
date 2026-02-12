core.set <- function(data.set, preset = NULL, coverage = 99, delta = 0.001, 
                       coverage_filename = "Coverage.csv",
                       Temp_file = "Temp.csv",
                       Coreset = "Coreset.csv"){
    
    # Check for data.table for fast I/O
    use_dt <- requireNamespace("data.table", quietly = TRUE)

    # Optimization: Work with matrix for speed
    if (is.data.frame(data.set)){
        counts <- as.matrix(data.set)
    } else {
        counts <- data.set
    }
    
    z <- Sys.time()
    print(z)
    
    # Store original data for final output and reference
    # data.set might be large, so we reference it via indices
    # We need access to columns by original index
    
    coverage <- as.numeric(coverage)
    delta <- delta
    
    # Track current indices of counts columns relative to original data.set
    # Initially 1:ncol
    current_indices <- seq_len(ncol(counts))
    cnames <- colnames(counts) # Original column names
    
    # Pre-calculate variable numbers (unique non-NA values per marker/row)
    # This is constant for markers
    var.num <- apply(counts, 1, function(x){length(unique(x[!is.na(x)]))})
    mpe <- sum(var.num - 1)
    
    result <- character(0)
    result.idx <- integer(0) # Indices in original data.set
    
    coverage.table <- NULL
    
    # Optimized overlap score function
    # Returns vector of scores for the input row
    overlap.score <- function(x){
        x_clean <- x[!is.na(x)]
        if(length(x_clean) == 0) return(numeric(length(x)))
        
        tbl <- table(x_clean)
        # Fast lookup
        idx <- match(x, names(tbl))
        sc <- as.numeric(tbl[idx])
        sc[is.na(sc)] <- 0
        return(sc)
    }
    
    prenum <- 0
    # Handle Preset
    if (!is.null(preset)){
        # Match preset names to original columns
        # Assuming preset contains column names
        preset.idx <- match(preset, cnames)
        preset.idx <- preset.idx[!is.na(preset.idx)]
        
        if(length(preset.idx) > 0){
            result <- cnames[preset.idx]
            result.idx <- preset.idx
            prenum <- length(result)
            
            # Update counts: remove alleles covered by preset
            # We must find where these preset columns are in 'counts'
            # At start, counts has all columns, so indices match
            
            # Extract preset data
            coreset_preset <- counts[, preset.idx, drop=FALSE]
            
            # Vectorized update of counts
            # For each preset sample, remove its alleles from counts
            for (i in 1:ncol(coreset_preset)){
                col_vec <- coreset_preset[, i]
                # Compare column-wise
                # We want to set elements in counts to NA if they match col_vec
                # Using which() handles NAs correctly (NA == NA is NA, which is dropped)
                counts[which(counts == col_vec)] <- NA
            }
            
            # Remove preset columns from counts and current_indices
            # To avoid re-selecting them
            # We need to find their indices in current counts
            # At this point, current_indices is 1..N
            # So we remove indices corresponding to preset.idx
            
            # Note: Removing columns changes indices.
            # Safe way: keep boolean mask or remove by name match?
            # Or just remove from current_indices and counts subsetting
            
            rem_mask <- current_indices %in% preset.idx
            counts <- counts[, !rem_mask, drop=FALSE]
            current_indices <- current_indices[!rem_mask]
            
            # Calculate initial coverage
            # Re-fetch from original data (or coreset_preset)
            # data.set might not be available as matrix if we didn't save it
            # But we have 'coreset_preset'
            
            ide.num <- apply(coreset_preset, 1, function(x){length(unique(x[!is.na(x)]))})
            coverage1 <- mean(ide.num/var.num*100)
            coverage.table <- rbind(coverage.table, c(0, "preset", coverage1, coverage1))
        }
    }
    
    # Main Loop
    for (idx in 1:mpe){
        if(ncol(counts) == 0) break
        
        cat(paste0(idx, "-th iteration starts at "), format(Sys.time()), "\n")
        
        # 1. Count non-NAs per sample
        not.na.counts <- colSums(!is.na(counts))
        if(max(not.na.counts) == 0) break # No more data
        
        # 2. Candidates
        candidate <- which(not.na.counts == max(not.na.counts))
        
        # 3. Calculate Overlap Scores
        # step0: Rows=Samples, Cols=Markers (because apply transposes)
        # Only compute for candidates? No, apply iterates rows (Markers).
        # We can't easily skip samples in apply(counts, 1, ...) without subsetting counts first.
        # But we need row-wise frequencies based on ALL remaining samples?
        # Logic: Frequency of allele in the *remaining* population.
        # So yes, use all current 'counts' columns.
        
        step0 <- apply(counts, 1, overlap.score)
        
        # 4. Mean Overlap
        # step0 is (Samples x Markers)
        # We need rowMeans for candidate rows in step0
        
        if (length(candidate) == 1){
            vals <- step0[candidate, ]
            overlap_val <- mean(vals[vals > 0], na.rm=TRUE) 
            overlap <- overlap_val
            names(overlap) <- colnames(counts)[candidate]
        } else {
            step01 <- step0[candidate, , drop=FALSE]
            # Replace 0 with NA for mean calculation if we want mean of non-zero matches?
            # Original: mean(x[!is.na(x)]). My score 0 for NA.
            # If 0 implies "no match" or "NA", we should exclude it?
            # Original logic: dummy[na.idx] <- 0.
            # So NA became 0.
            # Then mean(x[!is.na(x)]) included 0s?
            # No, 'idx' was match index.
            # If NA in data, match is NA.
            # dummy[NA] <- 0.
            # So NA became 0.
            # `overlap` was `mean(step01[!is.na(step01)])`.
            # `step01` had 0s. `!is.na` is TRUE for 0.
            # So 0s were INCLUDED.
            # So `rowMeans(step01)` is correct (as 0s are valid scores).
            overlap <- rowMeans(step01, na.rm=TRUE)
        }
        
        # 5. Select Best
        select_in_candidate <- which(overlap == max(overlap))
        
        final.rel.idx <- 0
        
        if (length(select_in_candidate) == 1){
            final.rel.idx <- candidate[select_in_candidate]
        } else {
            # Tie breaking using min var.num
            minvals <- numeric(length(select_in_candidate))
            for(i in seq_along(select_in_candidate)){
                s_idx <- candidate[select_in_candidate[i]]
                non_na <- which(!is.na(counts[, s_idx]))
                if(length(non_na) > 0) minvals[i] <- min(var.num[non_na]) else minvals[i] <- Inf
            }
            best_ties <- which(minvals == min(minvals))
            if(length(best_ties) == 1){
                final.rel.idx <- candidate[select_in_candidate[best_ties]]
            } else {
                final.rel.idx <- candidate[select_in_candidate[sample(best_ties, 1)]]
            }
        }
        
        # Map back to original index
        final.abs.idx <- current_indices[final.rel.idx]
        final.name <- cnames[final.abs.idx]
        
        result <- c(result, final.name)
        result.idx <- c(result.idx, final.abs.idx)
        
        cat("Selected:", final.name, "\n")
        
        # 6. Update Coverage
        # Use original data for coverage calculation
        # We need data.set. But we only have counts (which is modified).
        # We need to access original data for the selected sample.
        # We can keep a copy? Or read from file? 
        # Or, we can recover it? No.
        # We MUST keep original data or at least the selected columns.
        # But 'counts' has NAs.
        # Wait, we need to construct 'coreset' for output.
        # 'coreset' grows.
        # We can store 'coreset' separately and append columns.
        # But we need the *original* values, not the ones in 'counts' (which might have been NA-ed?).
        # Actually, 'counts' is modified by setting covered alleles to NA.
        # But the *selected* sample's column in 'counts' *at the moment of selection* 
        # might already have some NAs from previous iterations?
        # NO.
        # If a sample was not selected yet, its alleles are still there, 
        # UNLESS they were covered by OTHER samples?
        # The algorithm: `counts[idx1, i] = NA`.
        # This removes alleles from ALL samples if they match the selected one.
        # So `counts` represents "Uncovered Alleles".
        # But for the Output Coreset, we want the ORIGINAL alleles?
        # Usually yes.
        # Original code: `coreset <- data.frame(data.set[, result.idx])`.
        # It references `data.set`.
        # So we MUST keep `data.set` or a copy of it.
        # `counts <- data.set` (copy).
        # If `data.set` is large, keeping 2 copies is memory heavy.
        # But necessary if we want original values.
        # Can we rely on `counts`?
        # No, `counts` loses data.
        
        # So we assume `data.set` is available.
        # In my code `counts <- as.matrix(data.set)`. `data.set` is still there.
        # If `data.set` was passed as argument.
        
        # Construct current coreset for coverage calc
        if(is.matrix(data.set) || is.data.frame(data.set)){
             current_coreset <- data.set[, result.idx, drop=FALSE]
        } else {
             # Fallback if data.set was not kept (e.g. if we optimized memory)
             stop("Original data.set required for coverage calculation")
        }
        
        ide.num <- apply(current_coreset, 1, function(x){length(unique(x[!is.na(x)]))})
        coverage1 <- mean(ide.num/var.num*100)
        
        cat("Coverage is", coverage1, "%", "\n")
        
        # Diff calc
        if (idx == 1 & prenum == 0){
            dy <- coverage1
        } else {
            # coverage.table has previous rows
            prev_cov <- as.numeric(coverage.table[nrow(coverage.table), 3])
            dy <- coverage1 - prev_cov
        }
        cat("Difference is ", dy, "%", "\n")
        
        coverage.table <- rbind(coverage.table, c(idx, final.name, coverage1, dy))
        colnames(coverage.table) <- c("Iteration","Sample_name", "Coverage", "Difference")
        
        if (use_dt) {
            data.table::fwrite(as.data.frame(coverage.table), file = Temp_file)
        } else {
            write.csv(coverage.table, file = Temp_file, quote = FALSE)
        }

        if(coverage1 >= coverage | dy < delta){
            break
        }
        
        # 7. Update counts
        # Remove covered alleles
        # Use the column from 'counts' corresponding to selected sample?
        # No, 'counts' has NAs for covered alleles.
        # The selected sample covers the alleles it HAS.
        # Even if they were already covered?
        # If they were covered, they are NA in 'counts'.
        # So we only cover NEWly covered alleles?
        # Yes.
        # So we can use the current `counts[, final.rel.idx]`.
        
        selected_col <- counts[, final.rel.idx]
        
        # Vectorized Update
        counts[which(counts == selected_col)] <- NA
        
        # Remove the selected column from counts
        # to avoid selecting it again
        # We remove index `final.rel.idx`
        counts <- counts[, -final.rel.idx, drop=FALSE]
        current_indices <- current_indices[-final.rel.idx]
        
    }
    
    cat("GenoCore selects ", length(result), " element for core sets", "\n")
    cat("Running time is: \n")
    print(Sys.time() - z)
    
    coverage.table <- as.data.frame(coverage.table)
    
    # Final Output
    if(use_dt){
        data.table::fwrite(coverage.table, file = coverage_filename)
        data.table::fwrite(as.data.frame(data.set[, result.idx, drop=FALSE]), file = Coreset)
    } else {
        write.csv(coverage.table, file = coverage_filename, quote = FALSE, row.names=F)
        write.csv(data.set[, result.idx, drop=FALSE], file = Coreset, quote = FALSE, row.names = F)
    }
}
