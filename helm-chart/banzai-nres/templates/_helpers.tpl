{{/* vim: set filetype=mustache: */}}
{{/*
Expand the name of the chart.
*/}}
{{- define "banzai-nres.name" -}}
{{- default .Chart.Name .Values.nameOverride | trunc 63 | trimSuffix "-" -}}
{{- end -}}

{{/*
Create a default fully qualified app name.
We truncate at 63 chars because some Kubernetes name fields are limited to this (by the DNS naming spec).
If release name contains chart name it will be used as a full name.
*/}}
{{- define "banzai-nres.fullname" -}}
{{- if .Values.fullnameOverride -}}
{{- .Values.fullnameOverride | trunc 63 | trimSuffix "-" -}}
{{- else -}}
{{- $name := default .Chart.Name .Values.nameOverride -}}
{{- if contains $name .Release.Name -}}
{{- .Release.Name | trunc 63 | trimSuffix "-" -}}
{{- else -}}
{{- printf "%s-%s" .Release.Name $name | trunc 63 | trimSuffix "-" -}}
{{- end -}}
{{- end -}}
{{- end -}}

{{/*
Create chart name and version as used by the chart label.
*/}}
{{- define "banzai-nres.chart" -}}
{{- printf "%s-%s" .Chart.Name .Chart.Version | replace "+" "_" | trunc 63 | trimSuffix "-" -}}
{{- end -}}

{{/*
Common labels
*/}}
{{- define "banzai-nres.labels" -}}
app.kubernetes.io/name: {{ include "banzai-nres.name" . }}
helm.sh/chart: {{ include "banzai-nres.chart" . }}
app.kubernetes.io/instance: {{ .Release.Name }}
{{- if .Chart.AppVersion }}
app.kubernetes.io/version: {{ .Chart.AppVersion | quote }}
{{- end }}
app.kubernetes.io/managed-by: {{ .Release.Service }}
{{- end -}}


{{/*
Generate the PostgreSQL DB hostname
*/}}
{{- define "banzai-nres.dbhost" -}}
{{- if .Values.postgresql.fullnameOverride -}}
{{- .Values.postgresql.fullnameOverride | trunc 63 | trimSuffix "-" -}}
{{- else if .Values.useDockerizedDatabase -}}
{{- printf "%s-postgresql" .Release.Name -}}
{{- else -}}
{{- required "`postgresql.hostname` must be set when `useDockerizedDatabase` is `false`" .Values.postgresql.hostname -}}
{{- end -}}
{{- end -}}

{{- define "banzai-nres.rabbitmq" -}}
{{- if .Values.rabbitmq.fullnameOverride -}}
{{- .Values.rabbitmq.fullnameOverride | trunc 63 | trimSuffix "-" -}}
{{- else if .Values.useDockerizedRabbitMQ -}}
{{- printf "%s-rabbitmq" .Release.Name -}}
{{- else -}}
{{- required "`rabbitmq.hostname` must be set when `useDockerizedRabbitMQ` is `false`" .Values.rabbitmq.hostname -}}
{{- end -}}
{{- end -}}

{{/*
Define shared environment variables
*/}}
{{- define "banzai-nres.Env" -}}
- name: DB_HOST
  value: {{ include "banzai-nres.dbhost" . | quote }}
- name: RABBITMQ_HOST
  value: {{ include "banzai-nres.rabbitmq" . | quote }}
- name: DB_PASSWORD
  valueFrom:
    secretKeyRef:
      name: banzaiNresDBSecrets
      key: postgresPassword
- name: DB_USER
  value: {{ .Values.postgresql.postgresUsername }}
- name: DB_NAME
  value: {{ .Values.postgresql.postgresqlDatabase }}
- name: DB_ADDRESS
  value: psql://$(DB_USER):$(DB_PASSWORD)@$(DB_HOST)/$(DB_NAME)
- name: RABBITMQ_PASSWORD
  valueFrom:
    secretKeyRef:
      name: banzaiNresSecrets
      key: rabbitmq-password
- name: TASK_HOST
  value: rabbitmq://{{ .Values.rabbitmq.username }}@$(RABBITMQ_PASSWORD)$(RABBITMQ_HOST)/{{ .Values.rabbitmq.vhost }}
    secretKeyRef:
      name: banzaiNresSecrets
      key: TASK_HOST
- name: RETRY_DELAY
  value: 600000
- name: CALIBRATE_PROPOSAL_ID
  value: {{ .Values.CALIBRATE_PROPOSAL_ID }}
- name: OBSERVATION_PORTAL_URL
  value: {{ .Values.OBSERVATION_PORTAL_URL }}
- name: API_ROOT
  value:  {{ .Values.API_ROOT }}
- name: AUTH_TOKEN
  valueFrom:
    secretKeyRef:
      name: banzaiNresSecrets
      key: AUTH_TOKEN
- name: BUCKET
  value: {{ .Values.BUCKET }}
- name: AWS_ACCESS_KEY_ID
  valueFrom:
    secretKeyRef:
      name: banzaiNresSecrets
      key: AWS_ACCESS_KEY_ID
- name: AWS_SECRET_ACCESS_KEY
  valueFrom:
    secretKeyRef:
      name: banzaiNresSecrets
      key: AWS_SECRET_ACCESS_KEY
- name:
  value: {{ .Values.OPENTSDB_HOSTNAME }}
- name: BOSUN_HOSTNAME
  value: {{ .Values.BOSUN_HOSTNAME }}
- name: FITS_BROKER
  value: {{ .Values.FITS_BROKER }}
- name: FITS_EXCHANGE
  value: {{ .Values.FITS_EXCHANGE }}
- name: INGESTER_PROCESS_NAME
  value: {{ .Values.INGESTER_PROCESS_NAME }}
- name: POSTPROCESS_FILES
  value: "False"
- name: RAW_DATA_FRAME_URL
  value: {{ .Values.RAW_DATA_FRAME_URL }}
- name: RAW_DATA_AUTH_TOKEN
  valueFrom:
    secretKeyRef:
      name: banzaiNresSecrets
      key: RAW_DATA_AUTH_TOKEN
{{ if eq .Values.DO_METRICS 0 }}
- name: OPENTSDB_PYTHON_METRICS_TEST_MODE
  value: 1
{{- end -}}

{{- end -}}
