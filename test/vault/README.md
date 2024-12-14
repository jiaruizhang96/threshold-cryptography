
# README: How to Run a Vault Server with `n = 3` or `n = 5`

This guide provides detailed instructions for setting up and running a HashiCorp Vault server configured with either `n = 3` or `n = 5` unseal key shares. [Reference](https://developer.hashicorp.com/vault/tutorials/getting-started/getting-started-deploy).

## Prerequisites
1. **Vault Binary**: Ensure the `vault` binary is installed and accessible via your system's PATH.
2. **Configuration File**: Use the `config.hcl` file for Vault configuration.

## Steps to Set Up Vault

### 1. Stop Any Existing Vault Instance
If a Vault server is already running, terminate it to avoid conflicts.

```bash
pgrep -f vault | xargs kill
```

### 2. Clear Previous Data
Remove any existing history data to ensure a clean setup.

```bash
rm -r ./data
mkdir data
```

### 3. Start the Vault Server
Run the Vault server in the background using `nohup` and redirect logs to a file.

```bash
nohup sh -c "vault server -config=./config.hcl > ./log/vault.log 2>&1" > ./log/nohup.log &
```

### 4. Verify Vault Logs
Check the Vault server logs to ensure it started successfully.

```bash
cat ./log/vault.log
```

### 5. Export the Vault Address
Set the Vault address to match the server's API address.

```bash
export VAULT_ADDR=http://127.0.0.1:8200
```

### 6. Verify Vault Server Status
Check if the Vault server is running and ready for initialization.

```bash
vault status
```

## Initialize Vault

### 7. Initialize Vault with the Desired Key Shares and Threshold
Use `vault operator init` to initialize Vault with either `n = 3, t = 2` or `n = 5, t = 3`.

#### Example: `n = 3, t = 2`
```bash
vault operator init -key-shares=3 -key-threshold=2 > vault_init.txt
```

#### Example: `n = 5, t = 3`
```bash
vault operator init -key-shares=5 -key-threshold=3 > vault_init.txt
```

### 8. Save and Review Initialization Details
The output of the `vault operator init` command is saved to `vault_init.txt`. This file contains the unseal keys and the initial root token.

```bash
cat vault_init.txt
```

## Unseal Vault

### 9. Unseal Vault with Key Shares
Use the unseal keys from `vault_init.txt` to unseal Vault. You need to provide the number of shares specified by the threshold (`t`).

#### Example for `n = 3, t = 2`:
```bash
vault operator unseal $(grep 'Key 1:' vault_init.txt | awk '{print $NF}')
vault operator unseal $(grep 'Key 2:' vault_init.txt | awk '{print $NF}')
```

#### Example for `n = 5, t = 3`:
```bash
vault operator unseal $(grep 'Key 1:' vault_init.txt | awk '{print $NF}')
vault operator unseal $(grep 'Key 2:' vault_init.txt | awk '{print $NF}')
vault operator unseal $(grep 'Key 3:' vault_init.txt | awk '{print $NF}')
```

## Log in to Vault

### 10. Log in with the Initial Root Token
Use the root token from `vault_init.txt` to log in to Vault.

```bash
vault login $(grep 'Initial Root Token:' vault_init.txt | awk '{print $NF}')
```

## Enable Secrets Engine

### 11. Check Enabled Secrets Engines
List the currently enabled secrets engines.

```bash
vault secrets list
```

### 12. Enable the `kv` Secrets Engine
Enable the `kv` secrets engine at the path `secret/`.

```bash
vault secrets enable -path=secret kv
```

## Notes
- **Key Security**: Store the unseal keys and root token securely. 
- **Reinitialization**: If you need to reinitialize Vault, clear the data directory (`rm -r ./data`) and repeat these steps.
- **Logs**: Monitor `./log/vault.log` for any server-related issues.