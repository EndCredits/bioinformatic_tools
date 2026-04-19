{
  description = "TFBS Visualization Web Service";

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

      python = pkgs: pkgs.python313;
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
      pythonDeps = pkgs: with (python pkgs).pkgs; [ fastapi uvicorn python-multipart ];
    in
    {
      devShells = forAllSystems ({ system, pkgs }:
        {
          default = pkgs.mkShell {
            buildInputs = [
              (python pkgs)
              (rEnv pkgs)
            ] ++ (pythonDeps pkgs);

            shellHook = ''
              echo "TFBS Web Dev Shell"
              echo "  R:      $(R --version 2>&1 | head -1)"
              echo "  Python: $(python --version 2>&1)"
              echo ""
              echo "  Run:  R_SCRIPT=../visualize_tfbs.R uvicorn app:app --reload --port 6820"
            '';
          };
        });

      packages = forAllSystems ({ system, pkgs }:
        let
          app = pkgs.stdenv.mkDerivation {
            pname = "tfbs-web";
            version = "0.1.0";
            src = ./.;

            installPhase = ''
              mkdir -p $out/app/r $out/app/static
              cp $src/app.py $out/app/
              cp -r $src/static/* $out/app/static/
              cp ${../visualize_tfbs.R} $out/app/r/visualize_tfbs.R
            '';
          };
        in
        {
          # OCI image: nix build .#image && docker load < result
          #             or:  skopeo copy docker-archive:result docker://registry/tfbs-web:latest
          default = pkgs.dockerTools.buildLayeredImage {
            name = "tfbs-web";
            tag = "latest";
            contents = [
              (python pkgs)
              (rEnv pkgs)
              app
              pkgs.coreutils
              pkgs.bashInteractive
            ] ++ (pythonDeps pkgs);
            config = {
              Cmd = [
                "${(python pkgs).pkgs.uvicorn}/bin/uvicorn"
                "app:app"
                "--host" "0.0.0.0"
                "--port" "6820"
              ];
              Env = [ "R_SCRIPT=/app/r/visualize_tfbs.R" ];
              WorkingDir = "/app";
              ExposedPorts = { "6820/tcp" = {}; };
            };
          };
        });
    };
}
