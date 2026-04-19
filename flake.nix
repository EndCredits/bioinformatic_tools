{
  description = "Bioinformatic Tools";

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
    in
    {
      devShells = forAllSystems ({ system, pkgs }:
        {
          default = pkgs.mkShell {
            buildInputs = [
              pkgs.python313
              pkgs.R
              pkgs.git
            ];

            shellHook = ''
              echo "Bioinformatic Tools"
              echo "  Python: $(python --version 2>&1)"
              echo "  R:      $(R --version 2>&1 | head -1)"
            '';
          };
        });
    };
}
