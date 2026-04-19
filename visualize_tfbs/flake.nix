{
  description = "visualize_tfbs -- TFBS/Cis-Element Visualization";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" "x86_64-darwin" "aarch64-darwin" ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f {
        inherit system;
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = true;
        };
      });

      rPackages = pkgs: with pkgs.rPackages; [
        ggplot2
        dplyr
        tidyr
        viridis
        scales
        jsonlite
        yaml
      ];
      rEnv = pkgs: pkgs.rWrapper.override { packages = rPackages pkgs; };
    in
    {
      devShells = forAllSystems ({ system, pkgs }:
        {
          default = pkgs.mkShell {
            buildInputs = [
              (rEnv pkgs)
            ];

            shellHook = ''
              echo "visualize_tfbs Dev Shell"
              echo "  R: $(R --version 2>&1 | head -1)"
              echo ""
              echo "  Usage: Rscript visualize_tfbs.R [options] <fasta> <predictions> <output_prefix>"
              echo "         Rscript visualize_tfbs.R --plantcare [options] <tab_file> <fasta|-> <prefix>"
            '';
          };
        });

      packages = forAllSystems ({ system, pkgs }:
        {
          default = rEnv pkgs;
        });
    };
}
