/**
 * @license
 * Copyright 2025 Google LLC
 * SPDX-License-Identifier: Apache-2.0
 */

import { describe, it, expect, beforeEach, vi } from 'vitest';
import { H5adAnalyzerTool } from './h5ad-analyzer.js';
import * as fs from 'fs';

// Mock filesystem
vi.mock('fs');
vi.mock('child_process');

describe('H5adAnalyzerTool', () => {
  let tool: H5adAnalyzerTool;

  beforeEach(() => {
    tool = new H5adAnalyzerTool();
    vi.clearAllMocks();
  });

  describe('validateToolParams', () => {
    it('should reject empty h5adPath', () => {
      const result = tool.validateToolParams({
        h5adPath: '',
        analysisTask: 'test task',
      });
      expect(result).toBe('h5adPath is required');
    });

    it('should reject empty analysisTask', () => {
      const result = tool.validateToolParams({
        h5adPath: '/path/to/file.h5ad',
        analysisTask: '',
      });
      expect(result).toBe('analysisTask is required');
    });

    it('should reject non-h5ad files', () => {
      vi.mocked(fs.existsSync).mockReturnValue(true);
      const result = tool.validateToolParams({
        h5adPath: '/path/to/file.txt',
        analysisTask: 'test task',
      });
      expect(result).toBe('File must have .h5ad extension');
    });

    it('should reject non-existent files', () => {
      vi.mocked(fs.existsSync).mockReturnValue(false);
      const result = tool.validateToolParams({
        h5adPath: '/path/to/file.h5ad',
        analysisTask: 'test task',
      });
      expect(result).toContain('H5AD file not found');
    });

    it('should accept valid parameters', () => {
      vi.mocked(fs.existsSync).mockReturnValue(true);
      const result = tool.validateToolParams({
        h5adPath: '/path/to/file.h5ad',
        analysisTask: 'perform clustering',
      });
      expect(result).toBeNull();
    });
  });

  describe('getDescription', () => {
    it('should return formatted description', () => {
      const params = {
        h5adPath: 'test.h5ad',
        analysisTask: 'clustering analysis',
      };
      const result = tool.getDescription(params);
      expect(result).toBe(
        'Analyze h5ad file "test.h5ad" with task: "clustering analysis"',
      );
    });
  });

  describe('toolLocations', () => {
    it('should return h5ad file location', () => {
      const params = {
        h5adPath: 'test.h5ad',
        analysisTask: 'test',
      };
      const result = tool.toolLocations(params);
      expect(result).toHaveLength(1);
      expect(result[0].path).toContain('test.h5ad');
    });

    it('should include output directory if specified', () => {
      const params = {
        h5adPath: 'test.h5ad',
        analysisTask: 'test',
        outputDir: '/output',
      };
      const result = tool.toolLocations(params);
      expect(result).toHaveLength(2);
      expect(result[1].path).toContain('/output');
    });
  });

  describe('tool properties', () => {
    it('should have correct name and display name', () => {
      expect(tool.name).toBe('h5ad_analyzer');
      expect(tool.displayName).toBe('H5AD File Analyzer');
    });

    it('should have correct schema properties', () => {
      const schema = tool.schema;
      expect(schema.name).toBe('h5ad_analyzer');
      expect(schema.parameters?.properties).toHaveProperty('h5adPath');
      expect(schema.parameters?.properties).toHaveProperty('analysisTask');
      expect(schema.parameters?.required).toContain('h5adPath');
      expect(schema.parameters?.required).toContain('analysisTask');
    });
  });
});
